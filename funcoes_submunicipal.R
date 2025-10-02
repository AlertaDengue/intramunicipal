## Funções extras boletim submunicipal

# Função para normalizar nomes de bairros
normalizar_bairro <- function(nome) {
  nome %>%
    iconv(to = "ASCII//TRANSLIT") %>%
    toupper() %>%
    str_replace_all("[^A-Z0-9 ]", "") %>%
    str_squish()
}

# Função para ler dados de dengue (caso precise dos anos anteriores)
ler_dados_dengue <- function(caminho) {
  foreign::read.dbf(caminho) %>%
    filter(ID_MUNICIP == "420910", ID_MN_RESI == "420910") %>%  # Joinville
    select(NU_NOTIFIC, NU_ANO, SEM_NOT, SG_UF_NOT, ID_UNIDADE,
           ID_BAIRRO, NM_BAIRRO, ID_MN_RESI, CLASSI_FIN,
           DT_SIN_PRI, DT_DIGITA)
}

# Função para completar dados (bairros x semanas) - VERSÃO CORRIGIDA
completar_dados_bairros <- function(dados_bairros, bairros_unicos, se_var, nm_bairr_var, 
                                    group_var = "arbo", max_semana_epidemiologica = NULL) {
  
  # Load required libraries
  library(dplyr)
  library(tidyr)
  
  # Verify required columns exist
  required_cols <- c(se_var, nm_bairr_var, group_var, "notificações")
  missing_cols <- setdiff(required_cols, names(dados_bairros))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get unique arbovirus types
  arbo_unicos <- unique(dados_bairros[[group_var]])
  
  # Extract year from se_var (assuming format like 202401, 202402, etc.)
  dados_bairros <- dados_bairros %>%
    mutate(ano = as.numeric(substr(as.character(.data[[se_var]]), 1, 4)),
           semana = as.numeric(substr(as.character(.data[[se_var]]), 5, 6)))
  
  # Get all years present in data
  anos_unicos <- unique(dados_bairros$ano)
  
  # If max_semana_epidemiologica is provided, extract year and week from it
  if (!is.null(max_semana_epidemiologica)) {
    # max_semana_epidemiologica is in format YYYYWW (e.g., 202532)
    ano_atual <- as.numeric(substr(as.character(max_semana_epidemiologica), 1, 4))
    max_week_atual <- as.numeric(substr(as.character(max_semana_epidemiologica), 5, 6))
  } else {
    # Otherwise, determine current year and max week from data
    ano_atual <- max(anos_unicos)
    max_week_atual <- max(dados_bairros$semana[dados_bairros$ano == ano_atual], na.rm = TRUE)
  }
  
  # Create complete data for each year
  complete_data_list <- list()
  
  for (year in anos_unicos) {
    # Determine max week for this year
    if (year == ano_atual) {
      max_week <- max_week_atual
    } else {
      max_week <- 52  # Complete 52 weeks for previous years
    }
    
    # Create complete grid for this year
    year_grid <- expand_grid(
      !!nm_bairr_var := bairros_unicos,
      ano = year,
      semana = 1:max_week,
      !!group_var := arbo_unicos
    ) %>%
      mutate(!!se_var := as.numeric(paste0(ano, sprintf("%02d", semana))))
    
    complete_data_list[[as.character(year)]] <- year_grid
  }
  
  # Combine all years
  complete_data <- bind_rows(complete_data_list)
  
  # Join with original data and fill missing values with 0
  resultado <- complete_data %>%
    left_join(dados_bairros, by = c(nm_bairr_var, se_var, group_var, "ano", "semana")) %>%
    mutate(notificações = coalesce(notificações, 0)) %>%
    select(-ano, -semana) %>%  # Remove auxiliary columns
    arrange(!!sym(group_var), !!sym(nm_bairr_var), !!sym(se_var))
  
  return(resultado)
}

# Função para calcular alerta 
calcular_alerta <- function(id, nome, coluna_id, data, gtdist, meangt, sdgt) {
  bdi <- data %>%
    filter(.data[[coluna_id]] == id) %>%  # Filtra usando a coluna passada como string
    arrange(sem_not) %>%
    AlertTools::Rt(count = "casos", gtdist = gtdist, meangt = meangt, sdgt = sdgt)
  
  # Aplicar critérios de alerta
  crit.x.vector <- structure(as.character(params), names = names(params))
  criteriaU <- setCriteria(rule = params$codmodelo, values = crit.x.vector)
  y <- AlertTools::fouralert(bdi, crit = criteriaU)
  
  bdi$nivel <- y$indices$level
  bdi$receptivo <- y$indices$cytrue
  
  return(bdi)
}

# Função para classificar bairros por distrito
classificar_bairros_distrito <- function(sp_distritos, sp_bairros) {
  # Garantir mesmo CRS
  if(st_crs(sp_distritos) != st_crs(sp_bairros)) {
    sp_bairros <- st_transform(sp_bairros, st_crs(sp_distritos))
  }
  
  # Usar centroides para melhor precisão
  centroides <- st_centroid(sp_bairros)
  
  # Interseção espacial
  intersecao <- st_intersection(centroides, sp_distritos) %>%
    st_drop_geometry() %>%
    select(nome_bairr, id_bairro, distrito)
  
  # Completar com distritos mais próximos se necessário
  resultado <- sp_bairros %>%
    select(nome_bairr, id_bairro, geometry) %>%
    left_join(intersecao, by = c("nome_bairr", "id_bairro"))
  
  # Para bairros sem classificação, usar distrito mais próximo
  sem_distrito <- which(is.na(resultado$distrito))
  if(length(sem_distrito) > 0) {
    indices_proximos <- st_nearest_feature(
      st_centroid(resultado[sem_distrito,]), 
      sp_distritos
    )
    resultado$distrito[sem_distrito] <- sp_distritos$distrito[indices_proximos]
  }
  
  return(resultado %>% #select(-nome_bairr) %>%
           rename(geometry_distritos = geometry))
}

# Função otimizada para gerar mapas
generate_maps <- function(data, weeks_2_plotdistrito_nome, lim_range = NULL,
                          value_col = "p_inc100k", arbo = arbo, geom_col = "geometry",
                          shape = NULL, area = NULL) {
  
  # Calcular limites automaticamente se não fornecidos
  if (is.null(lim_range)) {
    filtered_data <- data %>% filter(sem_not %in% weeks_2_plot & arbo == !!arbo)
    lim_range <- range(pull(filtered_data, !!sym(value_col)), na.rm = TRUE)
  }
  
  # Definir os padrões de linha desejados
  line_types <- c("dotted", "dashed", "solid")
  
  # Gerar mapas usando map() para maior eficiência
  plots_list <- map(weeks_2_plot, ~ {
    week_data <- data %>% filter(sem_not == .x & arbo == !!arbo)
    
    # Calcular porcentagem de não alocados (se aplicável)
    percent_na <- 100 - (df_mun_ok %>% filter(sem_not == .x & arbo == !!arbo) %>% nrow() / 
                           df_mun %>% filter(sem_not == .x & arbo == !!arbo) %>% nrow() * 100)
    
    
    # Criar plot base
    p <- ggplot() +
      geom_sf(data = week_data %>% st_as_sf(), aes(fill = !!sym(value_col), geometry = !!sym(geom_col))) +
      scale_fill_distiller(palette = "Blues", direction = 1, limits = lim_range)
    
    # Adicionar camada de distritos com os padrões de linha especificados
    if (!is.null(shape)) {
      if ("geom" %in% names(shape) && !("geometry" %in% names(shape))) {
        shape <- shape %>% rename(geometry = geom)
      }
      
      if (!is.null(area)) {
        # Verificar quantos níveis únicos existem na coluna area
        area_levels <- unique(pull(shape, !!sym(area)))
        n_levels <- length(area_levels)
        
        # Criar um mapeamento de linhas para cada nível
        line_scale <- setNames(
          line_types[1:min(n_levels, length(line_types))], 
          area_levels[1:min(n_levels, length(line_types))]
        )
        
        p <- p + 
          geom_sf(data = shape %>% st_as_sf(), 
                  aes(geometry = geometry, linetype = !!sym(area)), 
                  linewidth = 0.5, alpha = 0, color = "black") +
          scale_linetype_manual(values = line_scale)
      } else {
        # Sem cor específica (apenas contorno)
        p <- p + 
          geom_sf(data = shape %>% st_as_sf(), 
                  aes(geometry = geometry), 
                  linewidth = 1, alpha = 0, color = "black", linetype = "solid")
      }
    }
    
    # Aplicar tema e labels
    p <- p +
      theme_pubclean(base_size = 10) +
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(),
            legend.key.width = unit(1, "cm")) +  # Aumenta espaço para visualizar linhas
      labs(title = paste("Semana de Notificação:", .x),
           subtitle = paste("Porcentagem de casos perdidos:", round(percent_na, 2), "%"))
    
    return(p)
  })
  
  return(plots_list)
}


create_incidence_plot <- function(nowcast_data, observed_data, title_suffix = NULL, 
                                  distrito_nome = NULL, SE = NULL,
                                  ylim_max = NULL, include_zoom = FALSE, pop = pop,
                                  convert_to_incidence = F) {
  
  # Transformar em incidência se solicitado
  if (convert_to_incidence) {
    observed_data <- observed_data %>%
      mutate(total_Y = total_Y / pop * 1e5)
    
    nowcast_data <- nowcast_data %>%
      mutate(across(c(Median, LI, LS), ~ .x / pop * 1e5))
    
    y_label <- "Incidência (casos por 100.000 hab.)"
  } else {
    y_label <- "Incidência (casos por 100.000 hab.)"
  }
  
  # Determinar o título baseado nos parâmetros fornecidos
  if (!is.null(title_suffix)) {
    plot_title <- title_suffix
  } else if (!is.null(distrito_nome) && !is.null(SE)) {
    plot_title <- paste0("Curva de incidência - ", distrito_nome, " até SE", SE)
  } else if (!is.null(distrito_nome)) {
    plot_title <- paste0("Curva de incidência - ", distrito_nome)
  } else {
    plot_title <- "Curva de incidência"
  }
  
  # Ensure we're comparing the same weeks
  common_weeks <- intersect(nowcast_data$ano_epi, observed_data$ano_epi)
  
  # Filter data to common weeks only
  nowcast_filtered <- nowcast_data %>% filter(ano_epi %in% common_weeks)
  observed_filtered <- observed_data %>% filter(ano_epi %in% common_weeks)
  
  # Apply the condition row-wise
  for (i in seq_along(common_weeks)) {
    if (nowcast_filtered$Median[i] < observed_filtered$total_Y[i]) {
      nowcast_filtered$Median[i] <- observed_filtered$total_Y[i]
    }
  }
  
  # Replace the original nowcast_data with the modified version
  nowcast_data <- nowcast_data %>%
    left_join(nowcast_filtered %>% select(ano_epi, Median), by = "ano_epi") %>%
    mutate(Median = coalesce(Median.y, Median.x)) %>%
    select(-Median.x, -Median.y)
  
  # max_val <- max(observed_data %>% 
  #                  filter(ano_epi >= paste0(ano_atual-1,substr(weeks_2_plot[3],5,6))) %>% 
  #                  pull(total_Y), na.rm = T)
  # 
  # min_val <- min(observed_data %>% 
  #                  filter(ano_epi >= paste0(ano_atual-1,substr(weeks_2_plot[3],5,6))) %>% 
  #                  pull(total_Y), na.rm = T)
  # 
  
  
  # Criar o gráfico principal
  p_main <- ggplot() +
    geom_segment(
      data = observed_data %>% 
        filter(ano_epi >= paste0(ano_atual-1,substr(weeks_2_plot[3],5,6))),
      aes(x = factor(ano_epi), xend = factor(ano_epi), 
          y = 0, yend = total_Y, color = factor(nivel)),
      linewidth = 0.7
    ) +
    # Linha observada
    geom_line(
      data = observed_data %>% 
        filter(ano_epi >= paste0(ano_atual-1,substr(weeks_2_plot[3],5,6))),
      aes(x = factor(ano_epi), y = total_Y, group = 1),
      color = "steelblue", linewidth = 1
    ) +
    # Nowcasting
    geom_line(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), y = Median, color = "Nowcasting", group = type),
      linewidth = 0.5, linetype = 2
    ) +
    # Ribbon nowcasting
    geom_ribbon(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"),
      aes(x = factor(ano_epi), ymin = LI, ymax = LS, fill = "Nowcasting", group = type),
      alpha = 0.2
    ) +
    geom_hline(yintercept = 76.4, col = "red", show.legend = F) +
    # Escalas
    scale_color_manual(
      name = NULL,
      values = c("Nowcasting" = "red3", "1" = "green", "2" = "yellow", 
                 "3" = "orange", "4" = "red3"),
      labels = c("Nowcasting" = "Nowcasting", "1" = "Baixo Risco", 
                 "2" = "Receptivo", "3" = "Transmissão", "4" = "Alta atividade"),
      breaks = c("Nowcasting", "1", "2", "3", "4")
    ) +
    scale_fill_manual(values = c("Nowcasting" = "red3"), guide = "none") +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(
      breaks = levels(factor(observed_data$ano_epi))[
        seq(1, nlevels(factor(observed_data$ano_epi)), by = 8)
      ]
    ) +
    # scale_y_continuous(
    #   limits= c(min_val,max_val)
    # ) +
    labs(
      x = "Semanas Epidemiológicas",
      y = "Incidência (casos por 100.000 hab.)",
      title = title_suffix
    )
    
    # Adicionar ylim se especificado
    if (!is.null(ylim_max)) {
      p_main <- p_main + ylim(0, ylim_max)
    }
  
  # Se não for incluir o zoom, retorne apenas o gráfico principal
  if (!include_zoom) {
    return(p_main)
  }
  
  # Preparar dados para o zoom (últimas 10 semanas)
  ultimas_semanas <- tail(sort(unique(observed_data$ano_epi)), 5)
  zoom_data <- observed_data %>% 
    filter(ano_epi %in% ultimas_semanas)
  
  # Criar gráfico de zoom
  p_zoom <- ggplot(zoom_data) +
    geom_segment(
      aes(x = factor(ano_epi), xend = factor(ano_epi),
          y = 0, yend = total_Y, color = factor(nivel)),
      linewidth = 0.8
    ) +
    geom_line(
      aes(x = factor(ano_epi), y = total_Y, group = 1),
      color = "steelblue", linewidth = 1.2
    ) +
    # Nowcasting
    geom_line(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"  & ano_epi %in% ultimas_semanas),
      aes(x = factor(ano_epi), y = Median, color = "Nowcasting", group = type),
      linewidth = 0.5, linetype = 2
    ) +
    # Ribbon nowcasting
    geom_ribbon(
      data = nowcast_data %>% dplyr::filter(type == "Nowcasting"  & ano_epi %in% ultimas_semanas),
      aes(x = factor(ano_epi), ymin = LI, ymax = LS, fill = "Nowcasting", group = type),
      alpha = 0.2
    )+
    geom_hline(yintercept = 72, col = "red", show.legend = F) +
    # Escalas
    scale_color_manual(
      name = NULL,
      values = c("Nowcasting" = "red3", "1" = "green", "2" = "yellow", 
                 "3" = "orange", "4" = "red"),
      labels = c("Nowcasting" = "Nowcasting", "1" = "Baixo Risco", 
                 "2" = "Receptivo", "3" = "Transmissão", "4" = "Alta atividade"),
      breaks = c("Nowcasting", "1", "2", "3", "4")
    ) +
    scale_fill_manual(values = c("Nowcasting" = "red3"), guide = "none") +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = "Últimas 5 semanas")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      title = element_text(size = 6),
      plot.background = element_rect(fill = "white", color = "gray"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(
      breaks = levels(factor(zoom_data$ano_epi))[
        seq(1, nlevels(factor(zoom_data$ano_epi)), by = 1)
      ]
    ) 
  
  # Combinar os gráficos usando patchwork
  final_plot <- p_main + 
    inset_element(p_zoom, left = 0.6, bottom = 0.3, right = 1, top = 1)
  
  return(final_plot)
}

# Função para criar gráfico de receptividade
create_receptivity_plot <- function(api_data, weeks_limit) {
  # Calcular fator de escala
  # scale_factor <- max(api_data$receptivo, na.rm = TRUE) / max(api_data$tempmin, na.rm = TRUE)
  
  api_data %>% 
    filter(SE >= paste0(ano_atual-1,substr(weeks_2_plot[3],5,6))) %>% 
    ggplot() +
    geom_rect(aes(xmin = as.numeric(factor(SE)) - 0.5, 
                  xmax = as.numeric(factor(SE)) + 0.5, 
                  ymin = -Inf, 
                  ymax = Inf, 
                  fill = tempmin >= 18), alpha = 0.7) +
    geom_line(aes(x = factor(SE), y = tempmin, group = 1), 
              color = "orange2", linewidth = 1) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"), guide = "none") +

    geom_hline(aes(yintercept = 18, color = "Limiar Favorável (18°C)"), 
               linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("red")) +
    scale_x_discrete(breaks = levels(factor(api_data$SE))[seq(1, nlevels(factor(api_data$SE)), by = 8)]) +
    scale_y_continuous(
      name = "Temperatura Mínima (°C)") +
    labs(x = "Semanas Epidemiológicas", 
         title = "Receptividade climática e Temperatura Mínima por Semana", 
         color = "") +
    theme_minimal() +
    theme(
      axis.title.y.right = element_blank(),  # Remove right axis title
      axis.text.y.right = element_blank(),   # Remove right axis text
      axis.title.y.left = element_text(color = "black"),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Função para criar gráfico de Rt (Criar condição para usar facet com distritos)
create_rt_plot <- function(api_data, weeks_limit, title_suffix, facet_by = NULL,
                           lwr = NULL, upr = NULL, sem_not = NULL) {
  
  # Set default column names if not provided
  if (is.null(sem_not)) sem_not <- "sem_not"
  
  # Check if required columns exist (only sem_not and Rt are mandatory)
  required_cols <- c(sem_not, "Rt")
  missing_cols <- required_cols[!required_cols %in% names(api_data)]
  if (length(missing_cols) > 0) {
    stop(paste("Colunas obrigatórias não encontradas nos dados:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check if optional columns exist when provided
  if (!is.null(lwr) && !lwr %in% names(api_data)) {
    warning(paste("Coluna", lwr, "não encontrada nos dados. Ribbon será ignorado."))
    lwr <- NULL
  }
  if (!is.null(upr) && !upr %in% names(api_data)) {
    warning(paste("Coluna", upr, "não encontrada nos dados. Ribbon será ignorado."))
    upr <- NULL
  }
  
  p <- api_data %>% 
    filter(!!sym(sem_not) <= as.numeric(as.character(weeks_limit))) %>% 
    ggplot() +
    geom_line(aes(x = factor(!!sym(sem_not)), y = Rt, group = 1), linetype = 1) +
    geom_hline(yintercept = 1, linetype = 2, col = "red") +
    scale_x_discrete(breaks = levels(factor(as.character(api_data[[sem_not]])))[seq(1, nlevels(factor(as.character(api_data[[sem_not]]))), by = 8)]) +
    labs(x = "Semanas Epidemiológicas", y = "Rt", title = title_suffix) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add ribbon if both lwr and upr are provided and exist
  if (!is.null(lwr) && !is.null(upr)) {
    p <- p + geom_ribbon(
      aes(x = factor(!!sym(sem_not)), 
          ymin = !!sym(lwr), 
          ymax = !!sym(upr), 
          group = if(!is.null(facet_by) && facet_by %in% names(api_data)) !!sym(facet_by) else 1),
      alpha = 0.7, fill = "grey"
    )
  }
  
  # Add faceting if specified and column exists
  if (!is.null(facet_by) && facet_by %in% names(api_data)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), ncol = 1)
  } else if (!is.null(facet_by)) {
    warning(paste("Variável", facet_by, "não encontrada nos dados. Facet ignorado."))
  }
  
  return(p)
}

# Função para criar painel de mapas
create_map_panel <- function(plots_list, legend_title = "Incidência (100.000 pessoas)") {
  # Extrair legenda
  legend_plot <- plots_list[[length(plots_list) - 2]] + 
    theme(legend.position = "top") + 
    labs(fill = legend_title)
  
  grob_legend <- ggplotGrob(legend_plot)
  legend_grob <- gtable_filter(grob_legend, "guide-box")
  
  # Criar painel
  ggdraw() +
    draw_plot(plots_list[[length(plots_list) - 2]] + theme(legend.position = "none"), 
              x = 0.08, y = 0.1, width = 0.27) +
    draw_plot(plots_list[[length(plots_list) - 1]] + theme(legend.position = "none"), 
              x = 0.65, y = 0.1, width = 0.27) +
    draw_plot(plots_list[[length(plots_list)]] + theme(legend.position = "none"), 
              x = 0.4, y = -0.3, width = 0.27) +
    draw_grob(legend_grob, x = 0.3, y = 0.05,width = 0.4) 
}


# Função para nowcasting por distrito
processar_nowcast_distrito <- function(dados_distritos,distritos_nomes, K = NULL, Dmax = NULL) {
  # Preparar dados
  # dados_distrito <- nowcastdg_dis %>% 
  #   left_join(df_bairros %>%
  #               select(nm_bairro_ref, distrito) %>% 
  #               distinct, 
  #             by = "nm_bairro_ref") %>% 
  #   select(nu_notific,sem_not,dt_digita,dt_sin_pri,distrito) %>% 
  #   #na.omit() %>% 
  #   filter(distrito == distrito_nome,
  #          year(dt_sin_pri) %in% c(ano_atual-2, ano_atual-1, ano_atual)) %>%
  #   select(dt_sin_pri, dt_digita, distrito)
  
  # dados_distrito[dados_distrito == ""] <- NA
  # dados_distrito <- na.omit(dados_distrito)
  
  # Gerar nowcasting
  df_nowcast <- nowcasting_inla(
    dataset = dados_distritos,
    data.by.week = TRUE,
    date_onset = dt_sin_pri,
    date_report = dt_digita,
    K = K, Dmax = Dmax
  )
  
  # Processar dados nowcast
  df_nowcast$total <- df_nowcast$total %>%
    mutate(
      ano = year(dt_event),
      epiweek = epical::epi_week(dt_event)[[1]],
      ano_epi = case_when(
        epiweek == 1 & lag(epiweek, default = 52) == 52 ~ 
          paste0(ano_atual, str_pad(epiweek, width = 2, side = "left", pad = "0")),
        TRUE ~ paste0(ano, str_pad(epiweek, width = 2, side = "left", pad = "0"))
      ),
      type = case_when(
        ano_epi > as.numeric(as.character(weeks_2_plot[3])) ~ "Forecast",
        TRUE ~ "Nowcasting"
      )
    )
  
  # Dados observados por semana
  dados_by_week <- dados_distritos %>%
    mutate(
      dt_event = dt_sin_pri,
      epiweek = epical::epi_week(dt_event)[[1]],
      ano = year(dt_event),
      ano_epi = paste0(ano, str_pad(epiweek, width = 2, pad = "0"))
    ) %>%
    count(ano_epi, dt_event, name = "total_Y") %>%
    filter(total_Y != 0 & ano_epi >= paste0(ano_atual-1, "01")) %>%
    group_by(ano_epi) %>%
    summarise(total_Y = sum(total_Y), .groups = "drop") %>%
    left_join(
      df_distritos_dengon %>% 
        filter(distrito == distritos_nomes) %>%
        mutate(sem_not = as.character(sem_not)) %>%
        select(sem_not, nivel),
      by = c("ano_epi" = "sem_not")
    )
  
  return(list(nowcast = df_nowcast, observado = dados_by_week))
}

# Função simples: contar quantas das últimas 5 semanas tiveram nível = 1
contar_semanas_nivel1 <- function(data, semana_atual, bairro, janela) {
  # Converter bairro para character
  bairro_char <- as.character(bairro)
  
  # Filtrar dados do bairro específico
  dados_bairro <- data[data$nm_bairro_ref == bairro_char, ]
  
  # Se não há dados, retorna 0
  if (nrow(dados_bairro) == 0) {
    return(0)
  }
  
  # Ordenar por semana
  dados_bairro <- dados_bairro[order(dados_bairro$sem_not), ]
  
  # Encontrar posição da semana atual
  pos_atual <- match(semana_atual, dados_bairro$sem_not)
  
  # Se não encontrar a semana atual, retorna 0
  if (is.na(pos_atual)) {
    return(0)
  }
  
  # Definir as últimas 'janela' semanas (incluindo a atual)
  inicio_janela <- max(1, pos_atual - janela + 1)
  fim_janela <- pos_atual
  
  # Extrair os níveis das últimas semanas
  niveis_janela <- dados_bairro$nivel[inicio_janela:fim_janela]
  
  # Contar quantas semanas tiveram nível = 1
  semanas_com_nivel1 <- sum(!is.na(niveis_janela) & niveis_janela == 1)
  
  return(semanas_com_nivel1)
}

# Função para contar número de semanas com alerta 3 e 4
calcular_semanas_alerta_alto <- function(data, semana_atual, bairro, janela) {
  bairro_char <- as.character(bairro)
  dados_bairro <- data[data$nm_bairro_ref == bairro_char, ]
  
  if (nrow(dados_bairro) == 0) return(0)
  
  dados_bairro <- dados_bairro[order(dados_bairro$sem_not), ]
  pos_atual <- match(semana_atual, dados_bairro$sem_not)
  
  if (is.na(pos_atual)) return(0)
  
  inicio_janela <- max(1, pos_atual - janela + 1)
  niveis <- dados_bairro$nivel[inicio_janela:pos_atual]
  
  # Contar níveis 3 e 4 (alerta alto)
  sum(!is.na(niveis) & niveis >= 3)
}
