pacman::p_load(
  tidyverse,purrr,lubridate,
  ncdf4,sf,terra,raster,INLA,reshape2
  )

## shapes

shape <- geobr::read_state(code_state = "SC")
shape_distrito <- read_sf("/home/ayrton/Documents/infodengue/intramunicipal-20250805T105042Z-1-001/intramunicipal/downscale_inla-20250805T123803Z-1-001/downscale_inla/Joishapes/joiDistritos/joiDistr.shp")

## funçõs para downscale
source("/home/ayrton/Documents/infodengue/intramunicipal-20250805T105042Z-1-001/intramunicipal/downscale_inla-20250805T123803Z-1-001/downscale_inla/funcoes_extras_downscale.R")

## carregando reanálise e ajustando dado

semana_interesse <- 1
semana_formatada <- sprintf("%02d", semana_interesse)

# Reading and cropping reanalysis ERA data
netcdf_reanalysis <- list.files(
  "/home/ayrton/Documents/infodengue/intramunicipal-20250805T105042Z-1-001/intramunicipal/downscale_inla-20250805T123803Z-1-001/downscale_inla/reanalysis(epiweek)/",
  pattern = paste0(semana_formatada, "\\.nc$"),
  full.names = TRUE
)



reanalysis_rasters <- lapply(
  netcdf_reanalysis,
  process_netcdf
)

initial_year <- tail(as.numeric(substr(terra::time(reanalysis_rasters[[1]]),1,4)))

df_reanalysis <- lapply(
  reanalysis_rasters,
  convert_to_dataframe
)

df_reanalysis <- bind_rows(
  df_reanalysis,
  .id = "year") %>% 
  tibble() %>% 
  mutate(
    year = initial_year - 1 + as.numeric(year),
    epiweek = semana_interesse) 

## mesh modelo

coords <- cbind(df_reanalysis$x, df_reanalysis$y)
boundary <- inla.nonconvex.hull(st_coordinates(shape)[, 1:2])

mesh <- inla.mesh.2d(
  loc = coords,
  boundary = shape,
  max.edge = c(0.5,1),
  offset = c(0.5, 1),
  cutoff = 0.9
)

## dados modelo

spde <- inla.spde2.pcmatern(
  mesh,
  alpha = 2, 
  constr = TRUE,
  prior.range = c(50, 0.5), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01)
)

timesn <- length(unique(df_reanalysis$year))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn
)

group <- as.numeric(factor(df_reanalysis$year, levels = unique(df_reanalysis$year), labels = 1:length(unique(df_reanalysis$year))))
A <- inla.spde.make.A(mesh = mesh, loc = coords, group = group)

stk.e <- inla.stack(
  tag = "est",
  data = list(y = df_reanalysis$values),
  A = list(1, A),
  effects = list(
    data.frame(
      b0 = rep(1, nrow(df_reanalysis)),
      year = factor(df_reanalysis$year, levels = unique(df_reanalysis$year))
    ), 
    s = indexs)
)

bb <- st_bbox(shape)
x <- seq(bb$xmin - 1, bb$xmax + 1, length.out = 500)
y <- seq(bb$ymin - 1, bb$ymax + 1, length.out = 500)
dp <- as.matrix(expand.grid(x, y))

p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y"))
st_crs(p) <- st_crs(shape)
ind <- st_intersects(shape, p)
dp <- dp[ind[[1]], ]

final_dp = NULL
for(i in unique(group)){
  temp_dp <- cbind(dp, Var3 = i)
  final_dp <- rbind(final_dp, temp_dp)
}

coop <- final_dp[, 1:2]
groupp <- final_dp[, 3]
final_dp <- final_dp %>% 
  data.frame() %>% 
  mutate(
    year = Var3,
    year = factor(year,
                  levels = unique(group),
                  labels = unique(df_reanalysis$year)
    )
  )

Ap <- inla.spde.make.A(mesh = mesh, loc = coop, group = groupp)

stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(
    data.frame(
      b0 = rep(1, nrow(final_dp)),
      year = final_dp$year
    ),
    s = indexs)
)

stk.full <- inla.stack(stk.e, stk.p)
rprior <- list(theta = list(prior = "pccor1", param = c(0.5, 0.9)))

formula <- y ~ 0 + b0 + year + f(s,
                                 model = spde, group = s.group,
                                 control.group = list(model = "ar1", hyper = rprior)
)

## rodando modelo

res <- inla(
  formula,
  data = inla.stack.data(stk.full),
  control.predictor = list(
    compute = TRUE,
    A = inla.stack.A(stk.full)
  )
)

## ajustando previsões

index <- inla.stack.index(stack = stk.full, tag = "pred")$data
final_dp <- data.frame(final_dp) %>% 
  set_names(c("x", "y", "time", "year"))

final_dp$pred_mean <- res$summary.fitted.values[index, "mean"]
final_dp$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
final_dp$pred_ul <- res$summary.fitted.values[index, "0.975quant"]

df_interpolate <- melt(final_dp,
                       id.vars = c("x", "y", "time", "year"),
                       measure.vars = c("pred_mean", "pred_ll", "pred_ul")
)

## cortando para município

# Filtrando

df_corte <- df_interpolate %>% 
  filter(x > -50 & y >= -27) #%>% 
  #filter(year %in% 2020:2024)

pontos_sf <- st_as_sf(
  df_corte,
  coords = c("x", "y"),
  crs = st_crs(shape_distrito)
)

pontos_dentro <- st_intersection(pontos_sf, shape_distrito)

df_filtrado <- as.data.frame(pontos_dentro)
df_filtrado <- cbind(df_filtrado, st_coordinates(pontos_dentro)) %>% 
  dplyr::select(-geometry)

## dados finais
df_salvar <-  df_filtrado  %>%
  pivot_wider(names_from = "variable",
              values_from = "value") %>%
  group_by(year,Distrito) %>% 
  #summarise(mean = mean(value)) %>% ungroup()
  summarise(mean_pred = mean(pred_mean),
            mean_ll = mean(pred_ll),
            pred_ul = mean(pred_ul)) %>% 
  mutate(SE = paste0(year,semana_interesse)) %>% 
  ungroup() %>% 
  dplyr::select(-year)


assign(paste0('res_epiweek', semana_interesse),df_salvar )
#head(res_epiweek)

write.csv(res_epiweek,"/home/ayrton/Documents/infodengue/intramunicipal-20250805T105042Z-1-001/intramunicipal/downscale_inla-20250805T123803Z-1-001/downscale_inla/resumo_boletim/")
