process_netcdf <- function(file) {
  # Carregar o arquivo NetCDF
  netcdf_data <- terra::rast(file)
  
  if (!identical(crs(netcdf_data), crs(shape))) {
    # Reprojetar o raster para o CRS do shapefile
    netcdf_data <- terra::project(netcdf_data, crs(shape))
  }
  
  # Recortar os dados usando o shapefile
  netcdf_crop <- terra::crop(netcdf_data, shape)
  netcdf_mask <- terra::mask(netcdf_crop, shape)
  
  # Retornar o raster recortado e mascarado
  return(netcdf_mask)
}


convert_to_dataframe <- function(r){
  # r = raw_rasters[[1]]
  df <- as.data.frame(r, xy = TRUE, cells = TRUE) %>% 
    gather(var, values, -c(cell, x, y))
  
  col_names <- colnames(df)
  
  df <- df %>% 
    mutate(
      values = values #- 273.15
    )
  
  return(df)
}



inla_spde <- function(df, spde, amat, apred, resp, covariaveis, df_covariaveis_interpoladas, raster_map, shape, title_label){
  
  # df = df_reanalysis
  # spde = spde
  # amat = Amat
  # apred = Apred
  # resp = "value"
  # covariaveis = NULL
  # df_covariaveis_interpoladas = df_mesh
  # raster_map = r
  # shape = shape
  # title_label = "t2m"
  
  inla_effects <- list(
    i = 1:spde$n.spde, # the spatial effect
    Intercept = rep(1, nrow(df))
  )
  for(i in covariaveis){
    inla_effects[[i]] <- df %>% pull(i)
  }
  
  vector_of_ones <- as.list(rep(1, length(inla_effects[covariaveis]) + 1))
  
  # Treinamento
  dat_stack <- inla.stack(data = list(resp = df %>% pull(resp)), # the response variable
                          A = c(amat, vector_of_ones), # the projection matrix
                          effects = inla_effects
  )
  
  # Efeito fixo
  df_covariaveis <- data.frame(inla_effects[covariaveis]) 
  if(nrow(df_covariaveis) > 0){
    modmat <- expand.grid(obter_sequencia_min_max(df_covariaveis))
    
    inla_effects_fixed <- list(
      Intercept = rep(1, nrow(modmat))
    )
    for(i in covariaveis){
      inla_effects_fixed[[i]] <- modmat %>% pull(i)
    }
    
    pred_stack_fixef <- inla.stack(data = list(resp = NA),
                                   A = vector_of_ones,
                                   effects = inla_effects_fixed,
                                   tag = "prd_fixef")
  }
  
  # Predicao
  inla_effects_pred <- list(
    Intercept = rep(1, nrow(df_covariaveis_interpoladas))
  )
  for(i in covariaveis){
    inla_effects_pred[[i]] <- df_covariaveis_interpoladas %>% pull(i)
  }
  
  pred_stack_alleff <- inla.stack(data = list(resp = NA),
                                  A = vector_of_ones,
                                  effects = inla_effects_pred,
                                  tag = "prd_alleff")
  
  # put all the stacks together
  if(nrow(df_covariaveis) > 0){
    all_stack <- inla.stack(dat_stack, pred_stack_fixef, pred_stack_alleff)
  }else{
    all_stack <- inla.stack(dat_stack, pred_stack_alleff)
  }
  
  formula_inla = paste("resp ~ -1 + Intercept +", paste(covariaveis, collapse=" + ") ,"+ f(i, model = spde)", collapse=" + ")
  formula_inla <- gsub("\\+\\s*\\+", "+", formula_inla)
  formula_inla <- gsub("\\+\\s*$", "", formula_inla)
  formula_inla = as.formula(formula_inla)
  
  m_inla <- inla(formula_inla,
                 data = inla.stack.data(all_stack),
                 control.predictor = list(A = inla.stack.A(all_stack), compute = TRUE),
                 quantiles = NULL)
  
  result_inla <- summary(m_inla)
  tabela_efeitos_fixos <- result_inla$fixed %>% 
    data.frame() %>% 
    mutate(
      lower_ci = mean - 1.96 * sd,
      upper_ci = mean + 1.96 * sd
    ) %>% 
    relocate(lower_ci, .after = mean) %>% 
    relocate(upper_ci, .after = lower_ci) %>% 
    knitr::kable(caption = paste0("Tabela: Efeitos fixos do modelo INLA."),
                 row.names = T,
                 col.names = c("Média", "IC inf.", "IC sup.", "SD", "Moda", "KLD"),
                 align = "c",
                 digits = 2,
                 format = "html") %>%
    kableExtra::kable_classic(full_width = F, html_font = "Cambria") %>% 
    row_spec(0, bold = T)
  
  id_fixef <- inla.stack.index(all_stack, "prd_fixef")$data
  
  if(!is.null(id_fixef)){
    modmat$resp <- m_inla$summary.fitted.values[id_fixef, "mean"]
    modmat$sd <- m_inla$summary.fitted.values[id_fixef, "sd"]
  }else{
    modmat <- data.frame(
      resp = m_inla$summary.fitted.values$mean,
      sd = m_inla$summary.fitted.values$mean
    )
  }
  
  id_alleff <- inla.stack.index(all_stack, "prd_alleff")$data
  
  df_resultado <- covariaveis_df %>% 
    mutate(
      pred = m_inla$summary.fitted.values[id_alleff, "mean"],
      pred_sd = m_inla$summary.fitted.values[id_alleff, "sd"],
      pred_lower_ci = with(covariaveis_df, pred - 1.96 * pred_sd/sqrt(nrow(covariaveis_df))),
      pred_upper_ci = with(covariaveis_df, pred + 1.96 * pred_sd/sqrt(nrow(covariaveis_df)))
    ) %>% 
    relocate(all_of(resp), .before = "pred")
  
  coordinates(df_resultado) <- ~ x + y
  gridded(df_resultado) <- TRUE
  
  covariaveis_raster <- brick(df_resultado)
  covariaveis_raster <- crop(covariaveis_raster, extent(shape))
  covariaveis_raster <- mask(covariaveis_raster, shape)
  proj4string(covariaveis_raster) <- crs(shape)
  
  coords <- coordinates(covariaveis_raster) %>% data.frame()
  covariaveis_df <- cbind(
    covariaveis_raster %>% 
      as.data.frame(), 
    coords
  ) %>% 
    na.omit()
  
  g_map <- covariaveis_df %>% 
    rename(response = pred) %>% 
    gerar_mapa(response = pred) +
    labs(title = title_label)
  
  result <- list(
    dados = covariaveis_df,
    mapa = g_map,
    efeitos_fixos = tabela_efeitos_fixos
  )
  
  return(result)
  
}

gg_interpolate <- function(sf, df_interpolate, var, 
                           mid_point = NA, max_value = NA, min_value = NA,
                           facet_var = NULL){
  
  # sf = shape
  # df_interpolate = df_interpolate
  # var = "pred_mean"
  # mid_point = NA
  # max_value = NA
  # min_value = NA
  
  df_temp <- df_interpolate %>% 
    mutate(value = ifelse(value < 0, NA, value))
  
  if(is.na(mid_point)){
    mid_point = round(median(df_temp$value, na.rm = T), 0)
  }
  
  if(is.na(max_value)){
    max_value = round(max(df_temp$value, na.rm = T) + 1, 0)
  }
  
  if(is.na(min_value)){
    min_value = round(min(df_temp$value, na.rm = T) - 1, 0)
  }
  
  df_temp <- df_temp %>% 
    filter(variable == var)
  
  g_map <- sf %>% 
    ggplot() + 
    geom_sf() + 
    coord_sf(datum = NA) +
    geom_tile(data = df_temp,
              aes(x = x, y = y, fill = value)) +
    labs(x = "Longitude", y = "Latitude", fill = "t2m") +
    scale_fill_viridis_c(
      na.value = "#ffffff",
      option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_value, max_value)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.title=element_text(size = 12), 
      legend.text=element_text(size = 12),
      axis.text = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.title.position = "top",
      legend.position = "right",
      # legend.direction = "horizontal",
      # legend.key.width = unit(2.5, "cm"),
      strip.text = element_text(size = 10)
    )
  
  if (!is.null(facet_var)) {
    g_map <- g_map +
      facet_wrap(as.formula(paste("~", facet_var)))
  }
  
  # df_sf <- st_as_sf(df_temp %>% 
  #                     dplyr::select(-c(time, variable)) %>% 
  #                     mutate(year = paste0("year_", year)) %>% 
  #                     pivot_wider(names_from = year, values_from = value),
  #                   coords = c("x", "y"), crs = 4326)
  # r <- raster(df_sf, res = c(0.01, 0.01)) 
  # raster_obj <- rasterize(df_sf, r)
  # plot(raster_obj)
  # 
  
  return(g_map)
  
}



