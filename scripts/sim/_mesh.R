# =====================================================================================================================
# STEP 1 NOTE: the true field is drawn on a mesh built from the real herd locations, and the analyst's mesh is built
# from the simulated herds, both with type "backbone" (Houaga's recipe). So, as in real life, the model's mesh is close
# to, but not exactly, the mesh the truth came from.
# =====================================================================================================================

# The SPDE mesh, built by ONE function for both the simulation and the fit: the fit mesh
# must equal the simulation mesh (note Part 4.10), so neither side may build it alone.
#   type "box"      — synthetic 80 x 80 km map (scripts/sim/00), herds not forced to be vertices
#   type "backbone" — Houaga's real-data recipe (1.Scripts_Spatial_Modelling.R, lines 157-198):
#                     inner non-convex hull 50 km around the herds, max.edge c(10, 40), cutoff 10.
#                     Houaga's outer boundary is the buffered Tanzania border; here a 200 km hull
#                     around the herds stands in for it (6 ranges of padding against edge effects).
# The data folder says which type to use in mesh_type.txt (missing file = "box").

build_mesh <- function(loc_herd, type = c("box", "backbone")) {
  type <- match.arg(type)
  if (type == "box") {
    return(fm_mesh_2d_inla(loc.domain = loc_herd, max.edge = c(8, 25), offset = c(12, 35), cutoff = 2))
  }
  bnd_inner <- fm_nonconvex_hull_inla(loc_herd, convex = 50)
  bnd_outer <- fm_nonconvex_hull_inla(loc_herd, convex = 200)
  fm_mesh_2d_inla(boundary = list(bnd_inner, bnd_outer), max.edge = c(10, 40), cutoff = 10)
}

mesh_type_of <- function(data_dir) {
  f <- file.path(data_dir, "mesh_type.txt")
  if (file.exists(f)) trimws(readLines(f, n = 1)) else "box"
}
