#!/usr/bin/env Rscript
# Usage: Rscript recife/fix_bongi_geometry.R SOURCE.RData OUTPUT.RData
# SOURCE must be the uncorrected upstream reference. OUTPUT may equal SOURCE.
suppressPackageStartupMessages(library(sf))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Expected SOURCE.RData and OUTPUT.RData')
source <- args[[1L]]
output <- args[[2L]]
original <- new.env(parent = emptyenv())
load(source, envir = original)
repaired <- new.env(parent = emptyenv())
for (name in ls(original, all.names = TRUE)) repaired[[name]] <- original[[name]]

stopifnot(identical(sort(ls(original)), sort(ls(repaired))))
stopifnot(nrow(original$rf_ref) == 94L,
          nrow(original$rf_rpa_shape) == 94L,
          nrow(original$rf_pop) == 94L,
          length(unique(original$rf_rpa_shape$rpa)) == 6L,
          sum(original$rf_pop$pop) == 1540874,
          !anyDuplicated(original$rf_ref$id_bairro),
          !anyDuplicated(original$rf_pop$id_bairro),
          !anyDuplicated(toupper(trimws(original$rf_ref$nome_bairr))),
          !anyDuplicated(toupper(trimws(original$rf_rpa_shape$bairro_nom))),
          identical(st_crs(original$rf_ref)$epsg, 4326L))

id <- 809L
ref_index <- which(original$rf_ref$id_bairro == id)
pop_index <- which(original$rf_pop$id_bairro == id)
stopifnot(length(ref_index) == 1L, length(pop_index) == 1L,
          original$rf_ref$nome_bairr[ref_index] == 'BONGI',
          original$rf_pop$bairro_nom[pop_index] == 'BONGI')
rpa_index <- which(original$rf_rpa_shape$bairro_nom == original$rf_pop$bairro_nom[pop_index])
stopifnot(length(rpa_index) == 1L)

canonical_fixed <- NULL
for (name in c('rf_ref', 'rf_rpa_shape')) {
  x <- original[[name]]
  index <- if (name == 'rf_ref') ref_index else rpa_index
  valid <- st_is_valid(x)
  stopifnot(identical(which(!valid), index), !any(st_is_empty(x)))
  before <- st_geometry(x)[index]
  if (name == 'rf_ref') {
    fixed <- st_make_valid(before)
  } else {
    stopifnot(identical(before[[1L]], st_geometry(original$rf_ref)[[ref_index]]))
    fixed <- canonical_fixed
  }
  stopifnot(length(fixed) == 1L, st_is_valid(fixed), !st_is_empty(fixed),
            as.character(st_geometry_type(fixed)) %in% c('POLYGON', 'MULTIPOLYGON'))
  # st_make_valid() returns a POLYGON here; keep the reference's MULTIPOLYGON type.
  fixed <- st_cast(fixed, 'MULTIPOLYGON')
  if (name == 'rf_ref') canonical_fixed <- fixed
  bbox_shift <- max(abs(as.numeric(st_bbox(before)) - as.numeric(st_bbox(fixed))))
  stopifnot(bbox_shift < 1e-7)
  before_area <- as.numeric(st_area(st_transform(st_set_crs(before, 4326), 31985)))
  after_area <- as.numeric(st_area(st_transform(st_set_crs(fixed, 4326), 31985)))
  relative_area_change <- abs(after_area - before_area) / before_area
  stopifnot(is.finite(relative_area_change), relative_area_change < 1e-4)
  st_geometry(x)[index] <- fixed
  repaired[[name]] <- x
  stopifnot(identical(st_crs(x), st_crs(original[[name]])),
            identical(st_drop_geometry(x), st_drop_geometry(original[[name]])),
            identical(class(x), class(original[[name]])),
            all(st_is_valid(x)), !any(st_is_empty(x)),
            all(vapply(setdiff(seq_len(nrow(x)), index),
                       function(i) identical(st_geometry(x)[[i]], st_geometry(original[[name]])[[i]]),
                       logical(1))))
  cat(sprintf('%s: Bongi %s -> %s; area %.3f -> %.3f m2 (change %.5f%%); max bbox shift %.9f degrees; all other geometries identical\n',
              name, as.character(st_geometry_type(before)),
              as.character(st_geometry_type(fixed)), before_area, after_area,
              100 * relative_area_change, bbox_shift))
}

for (name in setdiff(ls(original), c('rf_ref', 'rf_rpa_shape')))
  stopifnot(identical(original[[name]], repaired[[name]]))
stopifnot(identical(st_geometry(repaired$rf_ref)[[ref_index]],
                    st_geometry(repaired$rf_rpa_shape)[[rpa_index]]))
save(list = ls(repaired, all.names = TRUE), file = output, envir = repaired)

# Verify the serialized artifact in a fresh environment.
reloaded <- new.env(parent = emptyenv())
load(output, envir = reloaded)
stopifnot(identical(sort(ls(reloaded)), sort(ls(repaired))))
for (name in ls(repaired)) {
  equal <- if (name %in% c('fig1', 'fig2'))
    isTRUE(all.equal(reloaded[[name]], repaired[[name]]))
  else identical(reloaded[[name]], repaired[[name]])
  stopifnot(equal)
}
cat('Saved and reloaded all', length(ls(reloaded)), 'objects; 94 neighborhoods, 6 RPAs, population 1540874.\n')
