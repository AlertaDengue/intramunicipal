## pacotes

pacman::p_load(sf,sp,spdep,INLA,
               tidyverse,lubridate)

## espacial

rf_rpa <- rf_rpa_shape %>%
  group_by(rpa_nome) %>%
  sf::st_make_valid() %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()

rf_rpa_sp <- as(rf_rpa, "Spatial")
graph.file <- "RPA.graph"
nb2INLA(graph.file, poly2nb(rf_rpa_sp))

# 4. Visualizar a conectividade entre RPAs
plot(rf_rpa_sp, border = "gray")
plot(poly2nb(rf_rpa_sp), coordinates(rf_rpa_sp), col = "red", lwd = 2, add = TRUE)

## preparando dados

dados_individuais <- nowcastdg_rpa %>%
  select(
    rpa_nome,
    data_evento = dt_sin_pri,
    data_registro = dt_digita
  ) %>%
  # Filtrar dados válidos
  filter(!is.na(data_evento), !is.na(data_registro))

# Se QUISER agregar por semana (igual ao exemplo)
dados_semanais <- dados_individuais %>%
  mutate(
    # Criar semana epidemiológica (começando domingo)
    semana_evento = floor_date(data_evento, "week", week_start = 7),  # 7 = domingo
    semana_registro = floor_date(data_registro, "week", week_start = 7),
    # Delay em semanas
    delay_semanas = as.numeric(difftime(semana_registro, semana_evento, units = "weeks")),
    delay_semanas = round(delay_semanas)  # arredondar
  ) %>%
  # Agregar por distrito, semana e delay
  group_by(rpa_nome, semana_evento, delay = delay_semanas) %>%
  summarise(Y = n(), .groups = 'drop') %>%
  # Filtrar delays até 3-4 semanas (similar ao exemplo d0:d10 em semanas)
  filter(delay >= 0, delay <= 26) %>%
  # Completar combinações faltantes
  complete(rpa_nome, semana_evento, delay = 0:26, fill = list(Y = 0))

# Agora sim, seguir a mesma lógica do exemplo:
Today <- max(dados_semanais$semana_evento)  # 2 semanas de margem

# Filtrar até Today
dados_obs <- dados_semanais %>%
  filter(semana_evento <= Today)

# O mesmo loop para criar NAs (em semanas)
Dmax <- 10
for(d in Dmax:1) {
  dados_obs$Y[
    (dados_obs$semana_evento > (Today - 7*d)) & 
      (dados_obs$delay == d)
  ] <- NA
}

dados_obs %>% filter( is.na(Y) == T) %>% group_by(delay) %>% summarise(sum(is.na(Y)))

# Criar variável Time (semanas desde o início)
FirstDate <- min(dados_obs$semana_evento)
dados_obs <- dados_obs %>%
  mutate(
    Time = as.numeric(difftime(semana_evento, FirstDate, units = "weeks")) + 1
    )

dados_completos <- dados_semanais %>%  # antes do loop de NAs
  filter(semana_evento <= Today) %>%
  group_by(semana_evento) %>%
  summarise(Y_completo = sum(Y))

dados_observados <- dados_obs %>%  # após loop que criou NAs
  filter(semana_evento >= Today - 7*Dmax) %>%  # últimas Dmax semanas
  group_by(semana_evento) %>%
  summarise(Y_observado = sum(Y, na.rm = TRUE))

comparacao <- full_join(dados_completos, dados_observados, by = "semana_evento")

p0 <- ggplot(comparacao, aes(x = semana_evento)) +
  # Linha preta: dados completos (verdade)
  geom_line(aes(y = Y_completo, color = "Número real de casos")) +
  # Linha vermelha: dados observados (com atraso)
  geom_line(aes(y = Y_observado, color = "Casos notificados (com atraso)")) +
  labs(
    title = paste("Efeito do atraso de notificação - Today =", Today),
    x = "Semana do início dos sintomas (EpiWeekStart)",
    y = "Número de casos",
    color = "Legenda"
  ) +
  scale_color_manual(values = c("Número real de casos" = "black", 
                                "Casos notificados (com atraso)" = "red")) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.8, 0.9))

print(p0)

## modelo


## preparar espacial

# Criar mapeamento rpa_nome -> código numérico (1:n)
# IMPORTANTE: usar a mesma ordem do shapefile para consistência
rpa_ordem <- data.frame(
  rpa_nome = as.character(rf_rpa_sp@data$rpa_nome),
  GEOCOD = 1:nrow(rf_rpa_sp@data)  # ordem NATIVA do shapefile
)

# Juntar com dados_obs (garantir que nomes coincidem)
dados_obs <- dados_obs %>%
  mutate(rpa_nome = as.character(rpa_nome)) %>%
  left_join(rpa_ordem %>% select(rpa_nome, GEOCOD), by = "rpa_nome")

# Verificar se algum rpa_nome ficou sem GEOCOD
if(any(is.na(dados_obs$GEOCOD))) {
  stop("Existem rpa_nome em dados_obs que não estão no shapefile!")
}

# Criar índices para efeitos interativos
dados_obs <- dados_obs %>%
  mutate(
    # Delay.Region: interação delay x região (para efeitos iid)
    Delay.Region = paste(delay, GEOCOD, sep = "."),
    # Time.Region: interação tempo x região
    Time.Region = paste(Time, GEOCOD, sep = "."),
    # Para efeitos replicados no tempo (Time2 e Delay2)
    Time2 = Time,
    Delay2 = delay + 1,  # INLA precisa índices começando em 1
    # Índice para observações
    idx = 1:n()
  )

# Identificar índices com missing (para predição)
index.missing <- which(is.na(dados_obs$Y))

# Verificar estrutura
glimpse(dados_obs)

# criação do grafo

# O graph.file que você criou
graph.file <- "RPA.graph"

# Verificar o grafo
graph_info <- inla.read.graph(filename = graph.file)
print(paste("Número de regiões:", graph_info$n))
print(paste("Número de vizinhos:", sum(graph_info$nnbs)))
print("Matriz de adjacência (primeiras linhas):")
print(graph_info$graph[1:min(5, graph_info$n), 1:min(5, graph_info$n)])

####################################################
# 3. DEFINIR PRIORS (leo)
####################################################

half_normal_sd <- function(sigma) {
  return(
    paste0("expression:
              sigma = ", sigma, ";
              precision = exp(log_precision);
              logdens = -0.5*log(2*pi*sigma^2)-1.5*log_precision-1/(2*precision*sigma^2);
              log_jacobian = log_precision;
              return(logdens+log_jacobian);")
  )
}

####################################################
# 4. DEFINIR MODELOS COM ESTRUTURA ESPACIAL
####################################################

# Número máximo de delays
Dmax <- max(dados_obs$delay, na.rm = TRUE)

### Modelo 1: Efeitos principais + iid espacial + interação delay-região
modelo1 <- Y ~ 1 + 
  # Efeito temporal (RW1)
  f(Time, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) + 
  # Efeito do delay (RW1)
  f(delay, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(1)))  # prior mais fraco
  ) +
  # Efeito espacial iid (independente)
  f(GEOCOD, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  # Interação delay-região (iid)
  f(Delay.Region, model = "iid", constr = TRUE, 
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  )

### Modelo 2: Adicionar interação tempo-delay
modelo2 <- Y ~ 1 + 
  f(Time, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) + 
  f(delay, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(1)))
  ) +
  # Interação tempo-delay (efeito replicado)
  f(Time2, model = "rw1", constr = TRUE, 
    replicate = Delay2,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  f(GEOCOD, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  f(Delay.Region, model = "iid", constr = TRUE, 
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  )

### Modelo 3: Modelo completo com BYM espacial (mais sofisticado)
modelo3 <- Y ~ 1 + 
  f(Time, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) + 
  f(delay, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(1)))
  ) +
  f(Time2, model = "rw1", constr = TRUE, 
    replicate = Delay2,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  # Efeito espacial BYM (estruturado + não estruturado)
  f(GEOCOD, model = "bym", graph = graph.file, constr = TRUE,
    hyper = list(
      prec.unstruct = list(prior = half_normal_sd(0.1)),
      prec.spatial = list(prior = half_normal_sd(0.1))
    )
  ) +
  f(Delay.Region, model = "iid", constr = TRUE, 
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  )

### Modelo 4: Modelo completo + interação tempo-região
modelo4 <- Y ~ 1 + 
  f(Time, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) + 
  f(delay, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(1)))
  ) +
  f(Time2, model = "rw1", constr = TRUE, 
    replicate = Delay2,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  f(GEOCOD, model = "bym", graph = graph.file, constr = TRUE,
    hyper = list(
      prec.unstruct = list(prior = half_normal_sd(0.1)),
      prec.spatial = list(prior = half_normal_sd(0.1))
    )
  ) +
  f(Delay.Region, model = "iid", constr = TRUE, 
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  ) +
  # Interação tempo-região
  f(Time.Region, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = half_normal_sd(0.1)))
  )

## rodar modelo

# Escolher o modelo (comece com modelo3 que é o mais próximo do autor)
modelo_escolhido <- modelo3

# Rodar INLA
output <- inla(
  formula = modelo_escolhido,
  family = "nbinomial",
  data = as.data.frame(dados_obs),
  num.threads = 1,
  control.predictor = list(link = 1, compute = TRUE),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
  control.family = list(
    hyper = list(theta = list(prior = "loggamma", param = c(1, 0.1)))
  ),
  # Controlar a otimização
  control.inla = list(strategy = "adaptive", int.strategy = "eb")
)

# Resumo do modelo
summary(output)

## extrair resultados

# Fixed effects
print("=== FIXED EFFECTS ===")
print(output$summary.fixed)

# Hyperparameters
print("=== HYPERPARAMETERS ===")
print(output$summary.hyperpar)

# DIC e WAIC
print(paste("DIC:", round(output$dic$dic, 2)))
print(paste("WAIC:", round(output$waic$waic, 2)))

## amostragem nowscating

n_samples <- 1000
cat("Amostrando", n_samples, "da posteriori...\n")
samples_list <- inla.posterior.sample(n = n_samples, output)

cat("Amostrando valores faltantes...\n")
vector.samples <- lapply(
  X = samples_list,
  FUN = function(x, idx = index.missing) {
    mu <- exp(x$latent[idx])
    size <- x$hyperpar[1]  # parâmetro da binomial negativa
    rnbinom(n = length(idx), mu = mu, size = size)
  }
)

ultimas_semanas <- max(dados_obs$semana_evento) - 7 * 10

tibble.samples <- lapply(
  X = vector.samples,
  FUN = function(x) {
    data.aux <- dados_obs
    data.aux$Y[index.missing] <- x
    
    data.aux %>%
      filter(semana_evento >= ultimas_semanas) %>%
      group_by(rpa_nome, semana_evento) %>%
      summarise(Y = sum(Y), .groups = 'drop')
  }
)

state.pred <- bind_rows(tibble.samples, .id = "amostra")

## resumo por região
nowcasting_results <- state.pred %>%
  group_by(rpa_nome, semana_evento) %>%
  summarise(
    Media = mean(Y),
    Mediana = median(Y),
    LI = quantile(Y, probs = 0.025),
    LS = quantile(Y, probs = 0.975),
    Q25 = quantile(Y, probs = 0.25),
    Q75 = quantile(Y, probs = 0.75),
    # Probabilidade de exceder limiar (ajuste conforme seu contexto)
    # Prob_excedencia = mean(Y > 10),
    .groups = 'drop'
  )

dados_completos <- dados_semanais %>%
  filter(semana_evento <= Today) %>%
  group_by(rpa_nome, semana_evento) %>%
  summarise(Y_completo = sum(Y), .groups = 'drop')

dados_obs_agreg <- dados_obs %>%
  filter(!is.na(Y), semana_evento >= "2025-01-01") %>%
  group_by(rpa_nome, semana_evento) %>%
  summarise(Y_obs = sum(Y), .groups = 'drop')

## visu por região
regioes <- unique(nowcasting_results$rpa_nome)

lista_plots <- list()


for(reg in regioes) {
  
  p <- ggplot() +
    # Banda de incerteza do nowcasting
    geom_ribbon(
      data = filter(nowcasting_results, rpa_nome == reg),
      aes(x = semana_evento, ymin = LI, ymax = LS),
      fill = "gray70", alpha = 0.5
    ) +
    # Linha do nowcasting (mediana)
    geom_line(
      data = filter(nowcasting_results, rpa_nome == reg),
      aes(x = semana_evento, y = Mediana, color = "Nowcasting"),
      size = 1
    ) +
    # Linha dos dados observados com atraso
    geom_line(
      data = filter(dados_obs_agreg, rpa_nome == reg),
      aes(x = semana_evento, y = Y_obs, color = "Observados (com atraso)"),
      size = 0.8, linetype = "dashed"
    ) +
    # Linha dos dados completos (verdade)
    geom_line(
      data = filter(dados_completos %>% filter(semana_evento >= "2025-01-01"),
                    rpa_nome == reg),
      aes(x = semana_evento, y = Y_completo, color = "Completos (verdade)"),
      size = 0.8, linetype = "dotted"
    ) +
    # Linha vertical marcando Today
    geom_vline(xintercept = as.numeric(Today), linetype = "solid", color = "red", alpha = 0.5) +
    labs(
      title = paste("Nowcasting -", reg),
      subtitle = paste("Today =", Today, "| Dmax =", Dmax),
      x = "Semana do evento",
      y = "Número de casos",
      color = "Legenda"
    ) +
    scale_color_manual(
      values = c(
        "Nowcasting" = "black",
        "Observados (com atraso)" = "red",
        "Completos (verdade)" = "darkgreen"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  lista_plots[[reg]] <- p
  
}

plot_grid(plotlist = lista_plots, ncol = 2)

#### efeitos

# espacial BYM

# Extrair efeitos espaciais do modelo BYM
if("GEOCOD" %in% names(output$summary.random)) {
  
  # BYM tem dois componentes: estruturado (espacial) e não estruturado
  # No output, os IDs 1:n são o componente não estruturado
  # Os IDs (n+1):(2n) são o componente estruturado (espacial)
  n_regioes <- nrow(rpa_ordem)
  
  efeitos_espaciais <- output$summary.random$GEOCOD %>%
    mutate(
      tipo = ifelse(ID <= n_regioes, "não estruturado", "estruturado"),
      ID_orig = ifelse(ID <= n_regioes, ID, ID - n_regioes)
    ) %>%
    left_join(rpa_ordem %>% mutate(ID_orig = GEOCOD), by = "ID_orig")
  
  # Componente estruturado (espacial) para mapear
  efeitos_estruturados <- efeitos_espaciais %>%
    filter(tipo == "estruturado") %>%
    select(rpa_nome, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  
  # Juntar com shapefile
  rf_rpa_com_efeitos <- rf_rpa %>% as.data.frame() %>% 
    left_join(efeitos_estruturados)
  
  mapa1 <- ggplot(rf_rpa_com_efeitos) +
    geom_sf(aes(geometry=geometry, 
                fill = mean, title = "Efeito espacial")) +
    scale_fill_continuous(palette = "Reds") + 
    theme_minimal()
  
  print(mapa1)
  
  # Mapa da incerteza (desvio padrão)
  mapa2 <- ggplot(rf_rpa_com_efeitos) +
    geom_sf(aes(geometry=geometry, 
                fill = sd, title = "Efeito espacial")) +
    scale_fill_continuous(palette = "Reds") + 
    theme_minimal()
  
  print(mapa2)
}

# delay

if("delay" %in% names(output$summary.random)) {
  pDelay <- output$summary.random$delay %>%
    ggplot(aes(x = ID, y = `0.5quant`, 
               ymin = `0.025quant`, ymax = `0.975quant`)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, fill = "blue") +
    labs(title = "Efeito do delay",
         x = "Delay (semanas)", y = "Efeito") +
    scale_x_continuous(breaks = 0:Dmax) +
    theme_bw()
  
  print(pDelay)
}

# temporal

if("Time" %in% names(output$summary.random)) {
  FirstDate <- min(dados_obs$semana_evento)
  
  pTime <- output$summary.random$Time %>%
    mutate(Date = FirstDate + (ID - 1) * 7) %>%
    filter(Date >= "2025-01-01") %>% 
    ggplot(aes(x = Date, y = `0.5quant`, 
               ymin = `0.025quant`, ymax = `0.975quant`)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, fill = "green") +
    labs(title = "Efeito temporal",
         x = "Semana", y = "Efeito") +
    theme_bw()
  
  print(pTime)
}

