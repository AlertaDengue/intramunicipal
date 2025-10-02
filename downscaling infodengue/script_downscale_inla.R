r_pkgs <- c("ncdf4", "terra", "raster", "geobr", 
            "sf", "INLA", "kableExtra", "reshape2", 
            "tidyverse")

if(!"pacman" %in% rownames(installed.packages())) install.packages("pacman")
pacman::p_load(char = r_pkgs)
rm(r_pkgs)

source("./downscale/funcoes_extras_downscale.R")

# Parameters
#initial_year = 2000
#shape_mun <- read_municipality() |> filter(code_muni == "3304557")
#shape <- read_state(code_state = 33, year = 2020)

shape <- read_health_region() %>% 
  filter(abbrev_state == "SC" & name_health_region == "Nordeste")

# # Reading and cropping reanalysis ERA data
# netcdf_raw <- list.files("previsões(ecmwf51)/ok/", pattern = "\\.nc$", full.names = TRUE)
# raw_rasters <- lapply(netcdf_raw, process_netcdf)
# 
# df_raw <- lapply(raw_rasters, convert_to_dataframe)
# df_raw <- bind_rows(df_raw, .id = "year") %>%
#   tibble() %>%
#   mutate(
#     year = as.numeric(year) + (2000 - 1)
#   ) %>%
#   # filter(grepl("leadtime=1", var)) %>%
#   drop_na() %>%
#   group_by(cell, year, x, y) %>%
#   reframe(
#     values_mean = mean(values, na.rm = T),
#     values_sd = sd(values, na.rm = T)
#   ) %>%
#   dplyr::select(-cell)
# 
# min_val <- min(df_raw$values_mean, na.rm = T) - 1
# max_val <- max(df_raw$values_mean, na.rm = T) + 1
# 
# df_raw %>%
#   filter(year == 2024) %>%
#   ggplot() +
#   geom_tile(aes(x = x, y = y, fill = values_mean)) +
#   geom_point(aes(x = x, y = y), color = "red", size = 2) +
#   geom_sf(data = shape, fill = NA, color = "black", size = 0.5) +
#   coord_sf(datum = NA) +
#   scale_fill_viridis_c(option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_val, max_val)) +
#   labs(
#     title = "ERA5 t2m (raw)",
#     x = "Longitude",
#     y = "Latitude",
#     fill = "t2m"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = 12, face = "bold"),
#     legend.title=element_text(size = 12),
#     legend.text=element_text(size = 12),
#     axis.text = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     legend.title.position = "top",
#     legend.position = "right",
#     # legend.direction = "horizontal",
#     # legend.key.width = unit(2.5, "cm"),
#     strip.text = element_text(size = 10)
#   ) +
#   facet_wrap(~ year)

# Converting SpatRaster to data.frame and combining then into a single data.frame
# raw_dataframes <- lapply(raw_rasters, convert_to_dataframe)
# df_raw <- bind_rows(raw_dataframes, .id = "raster_id") %>% 
#   tibble() %>% 
#   mutate(
#     raster_id = as.numeric(raster_id) + (initial_year - 1)
#   ) %>% 
#   rename(year = raster_id)

# Reading and cropping reanalysis ERA data
netcdf_reanalysis <- list.files("./reanalysis(epiweek)/", pattern = "\\036.nc$", full.names = TRUE)
reanalysis_rasters <- lapply(netcdf_reanalysis, process_netcdf)
initial_year <- as.numeric(substr(terra::time(reanalysis_rasters[[1]]),1,4))

# Converting SpatRaster to data.frame and combining then into a single data.frame
df_reanalysis <- lapply(reanalysis_rasters, convert_to_dataframe)
df_reanalysis <- bind_rows(df_reanalysis, .id = "year") %>% 
     tibble() %>% 
     mutate(
       year = initial_year - 1 + as.numeric(year)
)


 min_val <- min(df_reanalysis$values, na.rm = T)
 max_val <- max(df_reanalysis$values, na.rm = T)

 df_reanalysis %>%
   ggplot(aes(x = x, y = y, fill = values)) +
   geom_raster() +
   geom_point(aes(x = x, y = y), color = "black", size = 0.5)+
   scale_fill_viridis_c(option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_val, max_val)) +
   theme_minimal() +
   labs(
     title = "ERA5 t2m (reanalysis)",
     x = "Longitude",
     y = "Latitude",
     fill = "t2m"
   ) +
   theme(legend.position = "right") +
   facet_wrap(~ year)

coords <- cbind(df_reanalysis$x, df_reanalysis$y)
boundary <- inla.nonconvex.hull(st_coordinates(shape_mun)[, 1:2])

mesh <- inla.mesh.2d(loc = coords,
                     boundary = shape,
                     max.edge = c(0.5,1),
                     offset = c(0.5, 1),
                     cutoff = 0.9)
# summary(mesh)
plot(mesh)
plot(st_geometry(shape), add = TRUE, border = "red")
summary(mesh)

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

res <- inla(formula,
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.full)
            )
)

result_inla <- summary(res)

result_inla$fixed %>% 
  data.frame() %>% 
  mutate(
    lower_ci = mean - 1.96 * sd,
    upper_ci = mean + 1.96 * sd
  ) %>% 
  relocate(lower_ci, .after = mean) %>% 
  relocate(upper_ci, .after = lower_ci) %>% 
  dplyr::select(mean:sd, mode, kld) %>% 
  knitr::kable(caption = paste0("Tabela: Efeitos fixos do modelo INLA."),
               row.names = T,
               col.names = c("Média", "IC inf.", "IC sup.", "DP", "Moda", "KLD"),
               align = "c",
               digits = 2,
               format = "html") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria") %>% 
  row_spec(0, bold = T)

result_inla$hyperpar %>% 
  data.frame() %>% 
  knitr::kable(caption = paste0("Tabela: Efeitos fixos do modelo INLA."),
               row.names = T,
               col.names = c("Média", "DP", "Q(2.75)", "Mediana", "Q(97.5)", "Moda"),
               align = "c",
               digits = 2,
               format = "html") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria") %>% 
  row_spec(0, bold = T)

list_marginals <- list(
  "b0" = res$marginals.fixed$b0,
  "2000" = res$marginals.fixed$year2000,
  "2001" = res$marginals.fixed$year2001,
  "2002" = res$marginals.fixed$year2002,
  "2003" = res$marginals.fixed$year2003,
  "2004" = res$marginals.fixed$year2004,
  "2005" = res$marginals.fixed$year2005,
  "2006" = res$marginals.fixed$year2006,
  "2007" = res$marginals.fixed$year2007,
  "2008" = res$marginals.fixed$year2008,
  "2009" = res$marginals.fixed$year2009,
  "2010" = res$marginals.fixed$year2010,
  "2011" = res$marginals.fixed$year2011,
  "2012" = res$marginals.fixed$year2012,
  "2013" = res$marginals.fixed$year2013,
  "2014" = res$marginals.fixed$year2014,
  "2015" = res$marginals.fixed$year2015,
  "2016" = res$marginals.fixed$year2016,
  "2017" = res$marginals.fixed$year2017,
  "2018" = res$marginals.fixed$year2018,
  "2019" = res$marginals.fixed$year2019,
  "2020" = res$marginals.fixed$year2020,
  "2021" = res$marginals.fixed$year2021,
  "2022" = res$marginals.fixed$year2022,
  "2023" = res$marginals.fixed$year2023,
  "precision Gaussian obs" =
    res$marginals.hyperpar$"Precision for the Gaussian observations",
  "range" = res$marginals.hyperpar$"Range for s",
  "stdev" = res$marginals.hyperpar$"Stdev for s",
  "rho" = res$marginals.hyperpar$"GroupRho for s"
)

marginals <- data.frame(do.call(rbind, list_marginals))
marginals$parameter <- factor(
  rep(names(list_marginals),
      times = sapply(list_marginals, function(x) if (is.null(dim(x))) length(x) else nrow(x))
  ),
  levels = names(list_marginals)
)

g_marginals <- marginals %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_line() +
  labs(x = "", y = "Density") + 
  theme_bw() +
  facet_wrap(~parameter, scales = "free", ncol = 3) 

# ggsave(plot = g_marginals,
#        filename = "figuras/g_marginals.png",
#        height = 20, width = 16, units = "cm", dpi = 300)

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


gg_mean <- gg_interpolate(sf = shape, df_interpolate = df_interpolate, var = "pred_mean", facet_var = "year") +
  labs(title = "ERJ - ERA5 (downscaled)",
       subtitle = "pred_mean")
gg_pred_ll <- gg_interpolate(sf = shape, df_interpolate = df_interpolate, var = "pred_ll", facet_var = "year")
gg_pred_ul <- gg_interpolate(sf = shape, df_interpolate = df_interpolate, var = "pred_ul", facet_var = "year")

gg_interpolate(sf = shape, df_interpolate = df_interpolate %>% filter(year == 2024), var = "pred_mean", facet_var = "year")

shape_mun <- read_municipality() %>% filter(code_muni == "4209102")

ggsave("./figures/ERJ_med09.png",gg_mean,height = 7,width = 10)


## Joinville municipalshapes

joiBairros <- read_sf("./Joishapes/joiBairros/joiBairros.shp") %>%
   select(id_bairro,geometry)
"~/Documents/vspinto/Lucas downscale/downscaling"
shape_mun <- read_sf("~/Documents/vspinto/Lucas downscale/downscaling/Joishapes/joiDistritos/joiDistr.shp")

ggplot() +
  geom_tile(data = df_interpolate %>% 
              filter(variable == "pred_mean"),
            aes(x = x, y = y, fill = value))+
  geom_sf(data = shape, aes(geometry = geom), alpha = 0, linewidth = 0.5) +
  #geom_sf(data = joiBairros, aes(geometry = geometry), alpha = 0, linewidth = 0.5) +
  #scale_fill_distiller(palette = "Reds")
  scale_fill_distiller(palette = "RdYlBu")  + 
  facet_wrap(~year)

# Converter df_interpolate para objeto sf (pontos)
pontos_sf <- st_as_sf(df_interpolate, coords = c("x", "y"), crs = st_crs(shape_mun))

# Realizar a interseção espacial
pontos_dentro <- st_intersection(pontos_sf, shape_mun)

# Converter de volta para dataframe se necessário
df_filtrado <- as.data.frame(pontos_dentro)
df_filtrado <- cbind(df_filtrado, st_coordinates(pontos_dentro)) %>% 
  select(-geometry)  # Remover a coluna geometry se não for mais necessária

# Visualizar o resultado

min_val <- min(df_filtrado$value, na.rm = T)
max_val <- max(df_filtrado$value, na.rm = T)

  ggplot(data = df_filtrado %>%
           filter(variable == "pred_mean"))+
  geom_tile(aes(x=X,y=Y, fill = value)) +
  geom_sf(data = shape_mun, aes(geometry = geometry), alpha = 0, linewidth = 0.5) +
    scale_fill_viridis_c(option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_val, max_val)) +
    facet_wrap(~year) +
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
      ) +
    labs(
          title = "RJ epiweek 01 - ERA5 t2m (downscaled)",
          subtitle = "pred_mean",
          x = "Longitude",
          y = "Latitude",
          fill = "°C"
        ) -> pp1; pp1

  ggplot(data = df_filtrado %>%
           filter(variable == "pred_ll"))+
    geom_tile(aes(x=X,y=Y, fill = value)) +
    geom_sf(data = shape_mun, aes(geometry = geometry), alpha = 0, linewidth = 0.5) +
    scale_fill_viridis_c(option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_val, max_val)) +
    facet_wrap(~year) +
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
    ) +
    labs(
      title = "Joinville - ERA5 t2m (downscaled)",
      subtitle = "Lower Interval",
      x = "Longitude",
      y = "Latitude",
      fill = "°C"
    ) -> pp2; pp2
  
  ggplot(data = df_filtrado %>%
           filter(variable == "pred_ul"))+
    geom_tile(aes(x=X,y=Y, fill = value)) +
    geom_sf(data = joiDistritos, aes(geometry = geometry), alpha = 0, linewidth = 0.5) +
    scale_fill_viridis_c(option = "turbo", breaks = scales::pretty_breaks(n = 5), limits = c(min_val, max_val)) +
    facet_wrap(~year) +
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
    ) +
    labs(
      title = " Joinville - ERA5 t2m (downscaled)",
      subtitle = "Upper Interval",
      x = "Longitude",
      y = "Latitude",
      fill = "°C"
    ) -> pp3; pp3

  ggsave("./figures/Join_med09.png",pp1,height = 7,width = 10)
  ggsave("./figures/Join_ll09.png",pp2,height = 7,width = 10)  
  ggsave("./figures/Join_ul09.png",pp3,height = 7,width = 10)  
  
  # 
df_filtrado  %>%
  pivot_wider(names_from = "variable",values_from = "value") %>% 
  group_by(year,Distrito) %>% 
  summarise(mean_pred = mean(pred_mean),
            mean_ll = mean(pred_ll),
            pred_ul = mean(pred_ul)) %>% 
  mutate(SE = paste0(year,"28")) %>% 
  ungroup() %>% 
  select(year) -> clim_downscaleJV

saveRDS(clim_downscaleJV,"~/Documents/vspinto/infodengue intramunicipal/infodengue_intramunicipal/Rt/climSE28.Rdata")  
