
# joga num diretório limpo e roda

### Pacotes
pacman::p_load(
  ecmwfr,tidyverse,
  purrr, ncdf4
)

### Área de interesse e chave (não funciona para municípios, necessita de downscaling pra isso)
wf_set_key(key = "b08f4374-52dc-4ff1-a8ff-7e67a8c2fa64")
shape <- geobr::read_health_region() %>% 
  filter(abbrev_state == "SC" & name_health_region == "Nordeste") # mudar aqui de acordo com o interesse
aoi_bb <- sf::st_bbox(shape)

# temperatura
### Request dos dados 
make_request <- function(year_month, variavel) {
  year <- year_month[1]
  month <- year_month[2]
  
  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = variavel,
    year = as.character(year),
    month = sprintf("%02d", as.numeric(month)),
    day = sprintf("%02d", 1:30), # todos os dias
    time = sprintf("%02d", 0:23), # todas as horas
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(aoi_bb[[2]], aoi_bb[[1]], aoi_bb[[4]], aoi_bb[[3]]),
    target = paste0("era5-",variavel, year, "_", month, "-demo.nc")
  )
  
  result <- wf_request(
    request = request,
    transfer = TRUE,
    path = ".",
    verbose = TRUE
  )
  
  # Pausa para evitar sobrecarregar a API
  Sys.sleep(2)
  
  return(result)
}

### Criar combinações de ano e mês
years <- 2014:2025  #especificar
months <- 7:9  # especificar
year_month_combinations <- expand.grid(year = years, month = months)

### Download temperatura
# Aplicar a função a todas as combinações
walk(
  array_branch(year_month_combinations, 1),
  make_request,
  variavel = "2m_temperature"
)

# a velocidade do vento no climate data store é divido de vários tipos, essa baixo é essa aqui:
# This parameter is the eastward component of the 10m wind. It is the horizontal speed of air moving
# towards the east, at a height of ten metres above the surface of the Earth, in metres per second.
# Care should be taken when comparing this parameter with observations, because wind observations
# vary on small space and time scales and are affected by the local terrain, vegetation and
# buildings that are represented only on average in the ECMWF Integrated Forecasting System (IFS).
# This parameter can be combined with the V component of 10m wind to give the speed and direction
# of the horizontal 10m wind.

# a vento v é a complementar que eles dizem precisar para estimar "the speed and direction
# of the horizontal 10m wind" Tem que ver como faz essa combinação, deve ter uma fórmula

### vento u ()
# Aplicar a função a todas as combinações
walk(
  array_branch(year_month_combinations, 1),
  make_request,
  variavel = "10m_u_component_of_wind"
)

### vento v
walk(
  array_branch(year_month_combinations, 1),
  make_request,
  variavel = "10m_v_component_of_wind"
)

### precipitação (pode baixar precipitação também)
# walk(
#   array_branch(year_month_combinations, 1),
#   make_request,
#   variavel = "total_precipitation"
# )

### Processando
### 1. função para calcular semana epi nos arquivos baixados
calcular_semana_epi <- function(data) {
  # Semana epidemiológica começa no domingo
  # Primeira semana do ano é a que contém 1º de janeiro
  
  # Trabalhar com cada data individualmente
  semanas <- numeric(length(data))
  
  for (i in 1:length(data)) {
    ano <- year(data[i])
    primeiro_jan <- as.Date(paste0(ano, "-01-01"))
    
    # Encontrar o primeiro domingo do ano ou anterior
    # wday() retorna 1 para domingo, 2 para segunda, etc.
    dia_semana_jan1 <- wday(primeiro_jan)
    primeiro_domingo <- primeiro_jan - (dia_semana_jan1 - 1)
    
    if (primeiro_domingo > primeiro_jan) {
      primeiro_domingo <- primeiro_domingo - 7
    }
    
    # Calcular número da semana
    dias_diff <- as.numeric(data[i] - primeiro_domingo)
    semanas[i] <- ceiling((dias_diff + 1) / 7)
  }
  
  return(semanas)
}
if(FALSE) {
### 2. função para quebrar em semana epidemiológica (atualmente gerando valores min)
processar_netcdf <- function(arquivo_entrada) {
  
  output_dir <- "reanalysis(epiweek)"
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Abrir arquivo de entrada
  nc_in <- nc_open(arquivo_entrada)
  
  # Ler dados
  var_names <- names(nc_in$var)
  arquivo <- var_names[[3]]
  long_name <- paste0("t$var$",arquivo,"$longname")
  temp_data <- ncvar_get(nc_in, arquivo)
  lon <- ncvar_get(nc_in, "longitude")
  lat <- ncvar_get(nc_in, "latitude")
  valid_time <- ncvar_get(nc_in, "valid_time")
  
  # Converter valid_time para datas
  # valid_time está em "seconds since 1970-01-01"
  datas <- as.POSIXct(valid_time, origin = "1970-01-01", tz = "UTC")
  datas_date <- as.Date(datas)
  
  # Calcular semanas epidemiológicas
  semanas_epi <- calcular_semana_epi(datas_date)
  
  # Agrupar por semana epidemiológica
  semanas_unicas <- unique(semanas_epi)
  
  cat("Processando", length(semanas_unicas), "semanas epidemiológicas\n")
  
  for (semana in semanas_unicas) {
    cat("Processando semana epidemiológica", semana, "\n")
    indices_semana <- which(semanas_epi == semana)
    
    if (arquivo == "t2m") {
      dados_agregados_semana <- apply(temp_data[, , indices_semana, drop = FALSE], c(1, 2), min, na.rm = TRUE)
      dados_agregados_semana <- dados_agregados_semana - 273.15
      unidade_saida <- "C"
    } else if (arquivo == "precip") {
      dados_agregados_semana <- apply(temp_data[, , indices_semana, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      dados_agregados_semana <- dados_agregados_semana * 1000
      unidade_saida <- "mm"
    } else if (arquivo == "u10") {
      dados_agregados_semana <- apply(temp_data[, , indices_semana, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      dados_agregados_semana <- dados_agregados_semana * 1000
      unidade_saida <- "ms"
    }  else if (arquivo == "v10") {
      dados_agregados_semana <- apply(temp_data[, , indices_semana, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      dados_agregados_semana <- dados_agregados_semana * 1000
      unidade_saida <- "ms"
    } 
    else {
      dados_agregados_semana <- apply(temp_data[, , indices_semana, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
      unidade_saida <- "desconhecida"
    }
  
  # # Criar arquivo para cada semana
  # for (semana in semanas_unicas) {
  #   cat("Processando semana epidemiológica", semana, "\n")
  #   
  #   # Indices of hours for this week
  #   indices_semana <- which(semanas_epi == semana)
  #   datas_semana <- datas_date[indices_semana]
  #       # # Calculate weekly statistic (average of daily values) with NA handling
    # t2m_semana <- apply(t2m_diario, c(1, 2), function(x) {
    #   if (all(is.na(x))) {
    #     return(NA)
    #   } else {
    #     return(mean(x, na.rm = TRUE))
    #   }
    # })
    # 
  #   # Get unique days in this week
  #   dias_unicos <- unique(datas_semana)
  #   n_dias <- length(dias_unicos)
  #   
  #   # Initialize array for daily statistics
  #   t2m_diario <- array(NA, dim = c(length(lon), length(lat), n_dias))
  #   
  #   # Calculate daily statistics for each day in the week
  #   for (d in 1:n_dias) {
  #     dia_atual <- dias_unicos[d]
  #     indices_dia <- indices_semana[datas_semana == dia_atual]
  #     
  #     if (length(indices_dia) > 0) {
  #       # MINIMUM: Daily minimum temperature with proper NA handling
  #       t2m_diario[, , d] <- apply(t2m[, , indices_dia, drop = FALSE], c(1, 2), function(x) {
  #         if (all(is.na(x))) {
  #           return(NA)
  #         } else {
  #           return(min(x, na.rm = TRUE))
  #         }
  #       })
        
        # MAXIMUM: Uncomment the lines below for daily maximum temperature
        # t2m_diario[, , d] <- apply(t2m[, , indices_dia, drop = FALSE], c(1, 2), function(x) {
        #   if (all(is.na(x))) {
        #     return(NA)
        #   } else {
        #     return(max(x, na.rm = TRUE))
        #   }
        # })
        
        # MEAN: Uncomment the lines below for daily mean temperature
        # t2m_diario[, , d] <- apply(t2m[, , indices_dia, drop = FALSE], c(1, 2), function(x) {
        #   if (all(is.na(x))) {
        #     return(NA)
        #   } else {
        #     return(mean(x, na.rm = TRUE))
        #   }
        # })
    #   }
    # }
    
    # # Calculate weekly statistic (average of daily values) with NA handling
    # t2m_semana <- apply(t2m_diario, c(1, 2), function(x) {
    #   if (all(is.na(x))) {
    #     return(NA)
    #   } else {
    #     return(mean(x, na.rm = TRUE))
    #   }
    # })
    # 
    
    # 
    # t2m_semana <- t2m_semana - 273.15
    
    # Calcular tempo representativo da semana (meio da semana)
    datas_semana <- datas_date[indices_semana]
    data_media <- mean(datas_semana)
    
    # Converter para "hours since 1900-01-01"
    referencia_1900 <- as.POSIXct("1900-01-01", tz = "UTC")
    horas_desde_1900 <- as.numeric(difftime(data_media, referencia_1900, units = "hours"))
    
    # Nome do arquivo de saída
    ano <- year(datas_semana[1])
    mes <- sprintf("%02d", month(datas_semana[1]))
    #nome_arquivo <- paste0("t2m_", ano, sprintf("%02d", semana), ".nc")
    output_file <- file.path(output_dir,paste0(arquivo,"_", ano, sprintf("%02d", semana), ".nc"))
    
    # Criar dimensões
    dim_lon <- ncdim_def("longitude", "degrees_east", lon, longname = "longitude")
    dim_lat <- ncdim_def("latitude", "degrees_north", lat, longname = "latitude")
    dim_time <- ncdim_def("time", "hours since 1900-01-01", horas_desde_1900, 
                          longname = "time", calendar = "gregorian")
    
    # Criar variável
    # var_t2m <- ncvar_def("t2m", "C", list(dim_lon, dim_lat, dim_time), 
    #                      missval = -32767, longname = "2 metre temperature")
    
    var_dados <- ncvar_def(arquivo, unidade_saida, list(dim_lon, dim_lat, dim_time),
                           missval = -32767, longname = long_name)
    
    nc_out <- nc_create(output_file, var_dados)
    ncvar_put(nc_out, var_dados, array(dados_agregados_semana, dim = c(dim(dados_agregados_semana), 1)))
    nc_close(nc_out)
    
    cat("Arquivo criado:", output_file, "\n")
  }
  
  # Fechar arquivo de entrada
  nc_close(nc_in)
  
  cat("Processamento concluído!\n")
}

era5_files <- list.files(pattern = "\\.nc")
lapply(era5_files, processar_netcdf)
}

# Função para processar o arquivo NetCDF
processar_netcdf <- function(arquivo_entrada) {
  
  output_dir <- "reanalysis(epiweek)"
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Abrir arquivo de entrada
  nc_in <- nc_open(arquivo_entrada)
  
  # Ler dados
  t2m <- ncvar_get(nc_in, "t2m")
  lon <- ncvar_get(nc_in, "longitude")
  lat <- ncvar_get(nc_in, "latitude")
  valid_time <- ncvar_get(nc_in, "valid_time")
  
  # Converter valid_time para datas
  # valid_time está em "seconds since 1970-01-01"
  datas <- as.POSIXct(valid_time, origin = "1970-01-01", tz = "UTC")
  datas_date <- as.Date(datas)
  
  # Calcular semanas epidemiológicas
  semanas_epi <- calcular_semana_epi(datas_date)
  
  # Agrupar por semana epidemiológica
  semanas_unicas <- unique(semanas_epi)
  
  cat("Processando", length(semanas_unicas), "semanas epidemiológicas\n")
  
  # Criar arquivo para cada semana
  for (semana in semanas_unicas) {
    cat("Processando semana epidemiológica", semana, "\n")
    
    # Indices das horas desta semana
    indices_semana <- which(semanas_epi == semana)
    
    # Calcular média da temperatura para esta semana
    #Média ao longo da dimensão temporal (3ª dimensão)
    #t2m_semana <- apply(t2m[, , indices_semana, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    
    t2m_semana <- apply(t2m[, , indices_semana, drop = FALSE], c(1, 2), function(x) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        return(min(x, na.rm = TRUE))
      }
    })
    
    t2m_semana <- t2m_semana - 273.15
    
    # Calcular tempo representativo da semana (meio da semana)
    datas_semana <- datas_date[indices_semana]
    data_media <- mean(datas_semana)
    
    # Converter para "hours since 1900-01-01"
    referencia_1900 <- as.POSIXct("1900-01-01", tz = "UTC")
    horas_desde_1900 <- as.numeric(difftime(data_media, referencia_1900, units = "hours"))
    
    # Nome do arquivo de saída
    ano <- year(datas_semana[1])
    mes <- sprintf("%02d", month(datas_semana[1]))
    #nome_arquivo <- paste0("t2m_", ano, sprintf("%02d", semana), ".nc")
    output_file <- file.path(output_dir,paste0("t2m_", ano, sprintf("%02d", semana), ".nc"))
    
    # Criar dimensões
    dim_lon <- ncdim_def("longitude", "degrees_east", lon, longname = "longitude")
    dim_lat <- ncdim_def("latitude", "degrees_north", lat, longname = "latitude")
    dim_time <- ncdim_def("time", "hours since 1900-01-01", horas_desde_1900, 
                          longname = "time", calendar = "gregorian")
    
    # Criar variável
    var_t2m <- ncvar_def("t2m", "C", list(dim_lon, dim_lat, dim_time), 
                         missval = -32767, longname = "2 metre temperature")
    
    # Criar arquivo NetCDF
    nc_out <- nc_create(output_file, var_t2m)
    
    # Escrever dados (adicionar dimensão temporal)
    ncvar_put(nc_out, var_t2m, array(t2m_semana, dim = c(dim(t2m_semana), 1)))
    # Fechar arquivo
    nc_close(nc_out)
    
    cat("Arquivo criado:", output_file, "\n")
  }
  
  # Fechar arquivo de entrada
  nc_close(nc_in)
  
  cat("Processamento concluído!\n")
}


era5_files <- list.files(pattern = "era5-")
lapply(era5_files, processar_netcdf)

