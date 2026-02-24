##### NOWCASTING ESPAÇO-TEMPORAL - VERSÃO GENERALIZADA

#========== pacotes

pacman::p_load(sf, sp, spdep, INLA, sn,
               tidyverse, lubridate, cowplot)

#========== PARÂMETROS DO USUÁRIO (edite aqui)

# Coluna que identifica a região no seu shapefile e nos dados individuais
col_regiao     <- "distrito"          # ex: "municipio", "distrito", "uf"

# Coluna da data do evento (início de sintomas)
col_data_evento    <- "dt_sin_pri"

# Coluna da data de registro/digitação
col_data_registro  <- "dt_digita"

# Shapefile já carregado como objeto sf
shape_sf       <- jv_distritos_shp      

# Dados individuais já carregados como data.frame
dados_raw      <- nowcastdg_distritos     

# Delay máximo (em semanas) a observar/prever
Dmax           <- 10

# Delay máximo a modelar (janela de delays armazenados)
Dmax_model     <- 26

# Nome do arquivo do grafo INLA
graph.file     <- "regiao.graph"

# Número de amostras da posteriori
n_samples      <- 1000

# Data de corte para visualização
data_inicio_viz <- as.Date("2025-01-01")

# Modelo a usar: 1, 2, 3 ou 4
modelo_escolhido_num <- 3

#========== FUNÇÕES AUXILIARES

half_normal_sd <- function(sigma) {
  paste0("expression:
              sigma = ", sigma, ";
              precision = exp(log_precision);
              logdens = -0.5*log(2*pi*sigma^2)-1.5*log_precision-1/(2*precision*sigma^2);
              log_jacobian = log_precision;
              return(logdens+log_jacobian);")
}

#========== PREPARAÇÃO ESPACIAL

# Dissolver polígonos por região e validar geometria
shape_regiao <- shape_sf %>%
  group_by(across(all_of(col_regiao))) %>%
  sf::st_make_valid() %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()

shape_sp <- as(shape_regiao, "Spatial")

# Criar grafo de vizinhança
nb2INLA(graph.file, poly2nb(shape_sp))

plot(shape_sp, border = "gray")
plot(poly2nb(shape_sp), coordinates(shape_sp),
     col = "red", lwd = 2, add = TRUE)

# Tabela de ordem nativa do shapefile -> índice GEOCOD
regiao_ordem <- data.frame(
  regiao = as.character(shape_sp@data[[col_regiao]]),
  GEOCOD = seq_len(nrow(shape_sp@data))
)
names(regiao_ordem)[1] <- col_regiao

#========== PREPARAÇÃO DOS DADOS

dados_individuais <- dados_raw %>%
  select(
    regiao      = all_of(col_regiao),
    data_evento = all_of(col_data_evento),
    data_registro = all_of(col_data_registro)
  ) %>%
  filter(!is.na(data_evento), !is.na(data_registro))

dados_semanais <- dados_individuais %>%
  mutate(
    semana_evento    = floor_date(data_evento,    "week", week_start = 7),
    semana_registro  = floor_date(data_registro,  "week", week_start = 7),
    delay_semanas    = round(as.numeric(difftime(semana_registro, semana_evento, units = "weeks")))
  ) %>%
  group_by(regiao, semana_evento, delay = delay_semanas) %>%
  summarise(Y = n(), .groups = "drop") %>%
  filter(delay >= 0, delay <= Dmax_model) %>%
  complete(regiao, semana_evento, delay = 0:Dmax_model, fill = list(Y = 0))

Today <- max(dados_semanais$semana_evento)

dados_obs <- dados_semanais %>%
  filter(semana_evento <= Today)

# Introduzir NAs para o período de nowcasting
for (d in Dmax:1) {
  dados_obs$Y[
    (dados_obs$semana_evento > (Today - 7 * d)) &
      (dados_obs$delay == d)
  ] <- NA
}

cat("NAs introduzidos por delay:\n")
print(dados_obs %>% filter(is.na(Y)) %>% group_by(delay) %>% summarise(n_NA = sum(is.na(Y))))

# Índice temporal (semanas desde a primeira data)
FirstDate <- min(dados_obs$semana_evento)
dados_obs <- dados_obs %>%
  mutate(Time = as.numeric(difftime(semana_evento, FirstDate, units = "weeks")) + 1)

# Juntar GEOCOD
dados_obs <- dados_obs %>%
  mutate(regiao = as.character(regiao)) %>%
  left_join(regiao_ordem %>% mutate(across(all_of(col_regiao), as.character)),
            by = setNames(col_regiao, "regiao"))

if (any(is.na(dados_obs$GEOCOD))) {
  stop("Há regiões em dados_obs que não constam no shapefile. Verifique a coluna '", col_regiao, "'.")
}

dados_obs <- dados_obs %>%
  mutate(
    Delay.Region = paste(delay, GEOCOD, sep = "."),
    Time.Region  = paste(Time,  GEOCOD, sep = "."),
    Time2        = Time,
    Delay2       = delay + 1,
    idx          = seq_len(n())
  )

index.missing <- which(is.na(dados_obs$Y))

#========== VISUALIZAÇÃO: EFEITO DO ATRASO

dados_completos_total <- dados_semanais %>%
  filter(semana_evento <= Today) %>%
  group_by(semana_evento) %>%
  summarise(Y_completo = sum(Y))

dados_observados_total <- dados_obs %>%
  filter(semana_evento >= Today - 7 * Dmax) %>%
  group_by(semana_evento) %>%
  summarise(Y_observado = sum(Y, na.rm = TRUE))

comparacao <- full_join(dados_completos_total, dados_observados_total, by = "semana_evento")

p0 <- ggplot(comparacao, aes(x = semana_evento)) +
  geom_line(aes(y = Y_completo,  color = "Número real de casos")) +
  geom_line(aes(y = Y_observado, color = "Casos notificados (com atraso)")) +
  labs(
    title = paste("Efeito do atraso de notificação - Today =", Today),
    x = "Semana do início dos sintomas", y = "Número de casos", color = NULL
  ) +
  scale_color_manual(values = c("Número real de casos" = "black",
                                "Casos notificados (com atraso)" = "red")) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.8, 0.9))

print(p0)

#========== GRAFO

graph_info <- inla.read.graph(filename = graph.file)
cat("Regiões no grafo:", graph_info$n, "\n")
cat("Total de vizinhanças:", sum(graph_info$nnbs), "\n")

#========== FÓRMULAS DOS MODELOS

f_modelo <- list(

  `1` = Y ~ 1 +
    f(Time,  model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(delay, model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(1)))) +
    f(GEOCOD, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(Delay.Region, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))),

  `2` = Y ~ 1 +
    f(Time,  model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(delay, model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(1)))) +
    f(Time2, model = "rw1", constr = TRUE, replicate = Delay2,
      hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(GEOCOD, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(Delay.Region, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))),

  `3` = Y ~ 1 +
    f(Time,  model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(delay, model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(1)))) +
    f(Time2, model = "rw1", constr = TRUE, replicate = Delay2,
      hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(GEOCOD, model = "bym", graph = graph.file, constr = TRUE,
      hyper = list(prec.unstruct = list(prior = half_normal_sd(0.1)),
                   prec.spatial  = list(prior = half_normal_sd(0.1)))) +
    f(Delay.Region, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))),

  `4` = Y ~ 1 +
    f(Time,  model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(delay, model = "rw1", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(1)))) +
    f(Time2, model = "rw1", constr = TRUE, replicate = Delay2,
      hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(GEOCOD, model = "bym", graph = graph.file, constr = TRUE,
      hyper = list(prec.unstruct = list(prior = half_normal_sd(0.1)),
                   prec.spatial  = list(prior = half_normal_sd(0.1)))) +
    f(Delay.Region, model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1)))) +
    f(Time.Region,  model = "iid", constr = TRUE, hyper = list(prec = list(prior = half_normal_sd(0.1))))
)

#========== RODAR INLA

output <- inla(
  formula  = f_modelo[[as.character(modelo_escolhido_num)]],
  family   = "nbinomial",
  data     = dados_obs,
  num.threads = 1,
  control.predictor = list(link = 1, compute = TRUE),
  control.compute   = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
  control.family    = list(hyper = list(theta = list(prior = "loggamma", param = c(1, 0.1)))),
  control.inla      = list(strategy = "adaptive", int.strategy = "eb")
)

#========== DIAGNÓSTICOS

cat("\n=== FIXED EFFECTS ===\n"); print(output$summary.fixed)
cat("\n=== HYPERPARAMETERS ===\n"); print(output$summary.hyperpar)
cat("DIC:", round(output$dic$dic, 2), "\n")
cat("WAIC:", round(output$waic$waic, 2), "\n")

#========== AMOSTRAS DA POSTERIORI

cat("Amostrando", n_samples, "da posteriori...\n")
samples_list <- inla.posterior.sample(n = n_samples, output)

cat("Imputando valores faltantes...\n")
vector.samples <- lapply(samples_list, function(x) {
  mu   <- exp(x$latent[index.missing])
  size <- x$hyperpar[1]
  rnbinom(n = length(index.missing), mu = mu, size = size)
})

ultimas_semanas <- Today - 7 * Dmax

tibble.samples <- lapply(vector.samples, function(x) {
  data.aux <- dados_obs
  data.aux$Y[index.missing] <- x
  data.aux %>%
    filter(semana_evento >= ultimas_semanas) %>%
    group_by(regiao, semana_evento) %>%
    summarise(Y = sum(Y), .groups = "drop")
})

state.pred <- bind_rows(tibble.samples, .id = "amostra")

nowcasting_results <- state.pred %>%
  group_by(regiao, semana_evento) %>%
  summarise(
    Media   = mean(Y),
    Mediana = median(Y),
    LI      = quantile(Y, 0.025),
    LS      = quantile(Y, 0.975),
    Q25     = quantile(Y, 0.25),
    Q75     = quantile(Y, 0.75),
    .groups = "drop"
  )

#========== DADOS AUXILIARES PARA PLOTS

dados_completos_reg <- dados_semanais %>%
  filter(semana_evento <= Today) %>%
  group_by(regiao, semana_evento) %>%
  summarise(Y_completo = sum(Y), .groups = "drop")

dados_obs_agreg <- dados_obs %>%
  filter(!is.na(Y), semana_evento >= data_inicio_viz) %>%
  group_by(regiao, semana_evento) %>%
  summarise(Y_obs = sum(Y), .groups = "drop")

#========== PLOTS POR REGIÃO

regioes <- unique(nowcasting_results$regiao)

lista_plots <- lapply(regioes, function(reg) {
  ggplot() +
    geom_ribbon(
      data = filter(nowcasting_results, regiao == reg),
      aes(x = semana_evento, ymin = LI, ymax = LS),
      fill = "gray70", alpha = 0.5
    ) +
    geom_line(
      data = filter(nowcasting_results, regiao == reg),
      aes(x = semana_evento, y = Mediana, color = "Nowcasting"), linewidth = 1
    ) +
    geom_line(
      data = filter(dados_obs_agreg, regiao == reg),
      aes(x = semana_evento, y = Y_obs, color = "Observados (com atraso)"),
      linewidth = 0.8, linetype = "dashed"
    ) +
    geom_line(
      data = filter(dados_completos_reg %>% filter(semana_evento >= data_inicio_viz), regiao == reg),
      aes(x = semana_evento, y = Y_completo, color = "Completos (verdade)"),
      linewidth = 0.8, linetype = "dotted"
    ) +
    geom_vline(xintercept = as.numeric(Today), linetype = "solid", color = "red", alpha = 0.5) +
    labs(
      title    = paste("Nowcasting -", reg),
      subtitle = paste("Today =", Today, "| Dmax =", Dmax),
      x = "Semana do evento", y = "Número de casos", color = NULL
    ) +
    scale_color_manual(values = c(
      "Nowcasting"               = "black",
      "Observados (com atraso)"  = "red",
      "Completos (verdade)"      = "darkgreen"
    )) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
})

names(lista_plots) <- regioes
plot_grid(plotlist = lista_plots, ncol = 2)

#========== EFEITOS DO MODELO

# Espacial (BYM ou iid)
if ("GEOCOD" %in% names(output$summary.random)) {
  n_regioes <- nrow(regiao_ordem)
  
  efeitos_espaciais <- output$summary.random$GEOCOD %>%
    mutate(
      tipo     = ifelse(ID <= n_regioes, "não estruturado", "estruturado"),
      ID_orig  = ifelse(ID <= n_regioes, ID, ID - n_regioes)
    ) %>%
    left_join(regiao_ordem %>% rename(ID_orig = GEOCOD), by = "ID_orig")
  
  efeitos_estruturados <- efeitos_espaciais %>%
    filter(tipo == "estruturado") %>%
    rename(regiao = all_of(col_regiao)) %>%
    select(regiao, mean, sd, `0.025quant`, `0.5quant`, `0.975qera pra ser dia 19, mas eu faltei e uant`)
  
  mapa_df <- shape_regiao %>%
    rename(regiao = all_of(col_regiao)) %>%
    left_join(efeitos_estruturados, by = "regiao")
  
  print(ggplot(mapa_df) +
    geom_sf(aes(geometry = geometry, fill = mean)) +
    scale_fill_viridis_c(option = "Reds", name = "Efeito médio") +
    labs(title = "Efeito espacial estruturado (BYM)") +
    theme_minimal())
  
  print(ggplot(mapa_df) +
    geom_sf(aes(geometry = geometry, fill = sd)) +
    scale_fill_viridis_c(option = "Reds", name = "Desvio padrão") +
    labs(title = "Incerteza do efeito espacial") +
    theme_minimal())
}

# Delay
if ("delay" %in% names(output$summary.random)) {
  print(
    output$summary.random$delay %>%
      ggplot(aes(x = ID, y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`)) +
      geom_line() + geom_ribbon(alpha = 0.5, fill = "blue") +
      scale_x_continuous(breaks = 0:Dmax_model) +
      labs(title = "Efeito do delay", x = "Delay (semanas)", y = "Efeito") +
      theme_bw()
  )
}

# Temporal
if ("Time" %in% names(output$summary.random)) {
  print(
    output$summary.random$Time %>%
      mutate(Date = FirstDate + (ID - 1) * 7) %>%
      filter(Date >= data_inicio_viz) %>%
      ggplot(aes(x = Date, y = `0.5quant`, ymin = `0.025quant`, ymax = `0.975quant`)) +
      geom_line() + geom_ribbon(alpha = 0.5, fill = "green") +
      labs(title = "Efeito temporal", x = "Semana", y = "Efeito") +
      theme_bw()
  )
}

