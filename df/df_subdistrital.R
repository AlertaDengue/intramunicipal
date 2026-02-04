
##esses df são criados com este script:
#df = dados até a semana mais recente
#df2 = dados linkados (acréscimo das colunas pop, Região de Saúde (RS) e Região Administrativa(RA))
#RS.summary = incidência por SE para cada RS
#RS.summary_complete = incidência por SE para cada RS (após verificar SE com 0 casos)
#nowcastDF.df = dados para nowcasting por RS




##packages
pacman::p_load(tidyverse,foreign,sf, cowplot,nowcaster, purr)

##load disease data and select few columns
df <- read.dbf("./data/DENGON_DF.dbf") %>%
      select(NU_NOTIFIC,NU_ANO,SEM_NOT,SEM_PRI,DT_NOTIFIC,DT_SIN_PRI,DT_DIGITA,
        SG_UF_NOT,SG_UF,ID_DISTRIT)


#load shapes and reference df for aggregate cases
load("./data/DF_subdistrital.RData")

##linking data
df2 <- df %>%
  left_join(regions %>% mutate(ID_DISTRITO  = as.factor(ID_DISTRITO)),
            by = c("ID_DISTRIT" = "ID_DISTRITO"))

head(df2)

df2 <-df2 %>%
  mutate(
    RS = ifelse(is.na(RS), "Não alocado", RS),
    RA = ifelse(is.na(RA), "Não alocado", RA)
)

missing <- df2 %>% filter(RS == "Não alocado") %>% nrow()
ok <- df2 %>% filter(RS != "Não alocado") %>% nrow()

missing_info <- print(pct_missing <- (missing/(ok+missing))*100)


# Calculating incidence ---------------------------------------------------


## RA
df2 %>% 
  group_by(RA, SEM_NOT) %>% 
  summarise(
    tot_casos = n(),
    tot_pop = sum(unique(pop)), 
    .groups = "drop"
  ) %>%
  mutate(inc = (tot_casos / tot_pop) * 1e5) -> RA.summary


### checking if every RA has the total number of epiweeks
### if not, include the missing epiweeks with 0 cases and incidence.
total_epiweeks <- as.numeric(substr(tail(df2$SEM_NOT,1),5,6))

missing_weeks <- RA.summary %>%
  group_by(RA) %>%
  summarise(num_weeks = n_distinct(SEM_NOT)) %>%
  filter(num_weeks < total_epiweeks)

### If RA with missing epiweeks
if (nrow(missing_weeks) > 0) {
  
first_epiweek <- as.numeric(as.character(head(df2$SEM_NOT, 1)))
last_epiweek <- as.numeric(as.character(tail(df2$SEM_NOT, 1)))
  

all_weeks <- data.frame(SEM_NOT = rep(first_epiweek:last_epiweek,
                                        each = length(unique(df2$RA))))

# fill missing epiweeks
RA.summary_complete <- RA.summary %>%
  group_by(RA) %>%
  complete(SEM_NOT = as.factor(first_epiweek:last_epiweek), fill = list(casos = 0, inc = 0)) %>%
  left_join(RAshape, by = "RA") %>% 
  ungroup()
}

##calculating incidence - RS
df2 %>% 
  group_by(RS, SEM_NOT) %>% 
  summarise(
    tot_casos = n(),
    tot_pop = sum(unique(pop)), 
    .groups = "drop"
 ) %>% 
  mutate(inc = (tot_casos / tot_pop) * 1e5) %>% 
  left_join(RSshape) -> RS.summary

# incidence maps ----------------------------------------------------------

weeks_2_plot <- tail(RA.summary_complete$SEM_NOT,3)  

## RA
limRA1 <- RA.summary_complete %>%
  filter(SEM_NOT %in% weeks_2_plot) %>% 
  pull(inc) %>% min()

limRA2 <- RA.summary_complete %>%
  filter(SEM_NOT %in% weeks_2_plot) %>% 
  pull(inc) %>% max()

plots_list_RA <- RA.summary_complete %>%
  filter(SEM_NOT %in% weeks_2_plot) %>%  
  split(.$SEM_NOT) %>%  
  map(~ ggplot(data = .x) +
        geom_sf(aes(fill = inc, geometry = geom)) +
        scale_fill_distiller(palette = "YlOrRd", direction = 1,
                             limits = c(limRA1, limRA2)) +
        ggpubr::theme_pubclean(base_size = 10) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank()) +
        labs(title = paste("Semana de Notificação:", unique(.x$SEM_NOT))))

### legend
legend <- plots_list_RA[[46]] + theme(legend.position = "top") + 
  labs(fill = "Incidência (100.000 pessoas)")
grob_p <- ggplotGrob(legend)
legend <- gtable::gtable_filter(grob_p, "guide-box")

### plot panel
ggdraw() +
  draw_plot(plots_list_RA[[46]] + theme(legend.position = "none"),x = 0.08, y =0.2,width = 0.4)+ 
  draw_plot(plots_list_RA[[47]] + theme(legend.position = "none"),x = 0.5, y =0.2,width = 0.4)+
  draw_plot(plots_list_RA[[48]] + theme(legend.position = "none"),x = 0.3, y =-0.24,width = 0.4) +
  draw_grob(legend, x = -0.02) -> RA_incidence_3epiweeks

## RS
limRS1 <- RS.summary %>%
  filter(SEM_NOT %in% weeks_2_plot) %>% 
  pull(inc) %>% min(na.rm = T)

limRS2 <- RS.summary %>%
  filter(SEM_NOT %in% weeks_2_plot) %>% 
  pull(inc) %>% max(na.rm = T)

plots_listRS <- RS.summary %>%
  filter(SEM_NOT %in% weeks_2_plot) %>% 
  split(.$SEM_NOT) %>%  
  map(~ ggplot(data = .x) +
        geom_sf(aes(fill = inc, geometry = geom)) +
        scale_fill_distiller(palette = "YlOrRd", direction = 1,
                             limits = c(limRS1, limRS2)) +
        ggpubr::theme_pubclean(base_size = 10) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank()) +
        labs(title = paste("Semana de Notificação:", unique(.x$SEM_NOT))))

### legend
legend <- plots_listRS[[46]] + theme(legend.position = "top") + 
  labs(fill = "Incidência (100.000 pessoas)")
grob_p <- ggplotGrob(legend)
legend <- gtable::gtable_filter(grob_p, "guide-box")

### plot panel
ggdraw() +
  draw_plot(plots_listRS[[46]] + theme(legend.position = "none"),x = 0.08, y =0.2,width = 0.4)+ 
  draw_plot(plots_listRS[[47]] + theme(legend.position = "none"),x = 0.5, y =0.2,width = 0.4)+
  draw_plot(plots_listRS[[48]] + theme(legend.position = "none"),x = 0.3, y =-0.24,width = 0.4) +
  draw_grob(legend, x = -0.02) -> RS_incidence_3epiweeks

# Nowcasting --------------------------------------------------------------
#data for nowcasting
nowcastDF.df <- df2 %>%  
  select(DT_DIGITA,DT_SIN_PRI,RS) %>% 
  arrange(RS) %>% 
  drop_na() %>% 
  filter(year(DT_SIN_PRI) %in% c(2023, 2024)) 

nowcastDF.df %>% mutate(RS_number = case_when(
  RS == "CENTRAL" ~ 1,
  RS == "CENTRO-SUL" ~ 2,
  RS == "Leste" ~ 3,
  RS == "NORTE" ~ 4,
  RS == "Não alocado" ~ 5,
  RS == "OESTE" ~ 6,
  RS == "SUDOESTE" ~ 7,
  RS == "SUL" ~ 8
)) -> nowcastDF.df



## loop for nowcasting
unique_RS <- unique(nowcastDF.df$RS)


generate_plot <- function(rs_value) {
  # Filtrar os dados para a RS específica
  filtered_data <- nowcastDF.df %>% filter(RS == rs_value)
  
  # Gerar o nowcasting
  df_nowcast <- nowcasting_inla(
    dataset = filtered_data,
    data.by.week = TRUE,
    date_onset = DT_SIN_PRI,
    date_report = DT_DIGITA,
    K = 3
  )
  
  # Ajustar os dados de nowcast
  df_nowcast$total <- df_nowcast$total %>% 
    mutate(epiweek = lubridate::epiweek(dt_event)) %>% 
    mutate(type = case_when(
      epiweek <= total_epiweeks ~ "Nowcasting",
      epiweek > total_epiweeks ~ "Forecast"
    ))
  
  # Dados observados agregados por semana
  dados_by_week <- df_nowcast$data %>%  
    ungroup() %>%  
    select(dt_event, Y) %>%  
    group_by(dt_event) %>%  
    summarize(total_Y = sum(Y, na.rm = TRUE), .groups = "drop") %>% 
    mutate(epiweek = epiweek(dt_event))
  
  # Criar o gráfico
  plot <- ggplot() + 
    geom_line(data = df_nowcast$total,
              aes(x = epiweek, y = Median, col = type)) +
    geom_ribbon(data = df_nowcast$total,
                aes(x = epiweek, ymin = LI, ymax = LS, col = type, fill = type), 
                alpha = 0.2, show.legend = F) +
    geom_line(data = dados_by_week %>% filter(epiweek <= total_epiweeks),
              aes(x = epiweek, y = total_Y, col = "Observed")) +
    scale_x_continuous(breaks = seq(min(dados_by_week$epiweek), max(df_nowcast$total$epiweek), by = 6)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45)) +
    scale_color_manual(values = c('red3', 'grey50', "black"), name = '') +
    scale_fill_manual(values = c('red3', 'grey50'), name = '') +
    labs(
      x = 'Semana epidemiológica', 
      y = 'Nº de Casos', 
      title = paste("Região de Saúde -", rs_value)
    )
  
  # Adicionar subtítulo condicional
  if (rs_value == "Não alocado") {
    plot <- plot + labs(subtitle = paste(round(missing_info,2),"% casos não alocados em nenhuma RS"))
  }
  
  return(plot)
}

# Gerar os gráficos para cada RS e salvar em uma lista
plots_list_nowcasting <- map(unique_RS, ~ generate_plot(.x))

ggpubr::ggarrange(plots_list_nowcasting[[1]],
          plots_list_nowcasting[[2]],
          plots_list_nowcasting[[3]],
          plots_list_nowcasting[[4]],
          plots_list_nowcasting[[5]],
          plots_list_nowcasting[[6]],
          plots_list_nowcasting[[7]],
          plots_list_nowcasting[[8]],
          ncol = 4, nrow = 2, common.legend = T) -> RS_nowcasting


# Results -----------------------------------------------------------------


RA_incidence_3epiweeks
RS_incidence_3epiweeks
RS_nowcasting
