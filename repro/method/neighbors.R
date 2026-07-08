# M0 sparse neighbor engine for lattice spatial-transcriptomics data.
# ST spots sit on an integer lattice, so the "Euclidean dist <= 1" neighbor rule
# used in R/dsgd.R and R/prediction.R reduces to the 4 rook neighbors and can be
# found by coordinate hashing in O(N*k) -- no kd-tree / extra deps needed.
suppressMessages(library(Matrix))

# Rook offsets: exactly the integer-lattice points with Euclidean distance <= 1.
.ROOK <- list(c(1L,0L), c(-1L,0L), c(0L,1L), c(0L,-1L))

# Precomputed neighbor sum  c_i = sum_{j in N(i)} y_j   (label-matched if given).
# During fitting the neighbor labels y are OBSERVED data, so c is a constant vector
# -- this is the whole point of M0: the dense NxN loop recomputes this constant.
neighbor_sum <- function(coords, y, label = NULL, radius = 1) {
  stopifnot(nrow(coords) == length(y))
  n  <- length(y)
  ix <- as.integer(round(coords[, 1]))
  iy <- as.integer(round(coords[, 2]))
  if (is.null(label)) label <- rep(1L, n)
  # hash "x,y,label" -> row index; look up each rook neighbor
  key <- paste(ix, iy, label, sep = ",")
  idx <- seq_len(n)
  lookup <- setNames(idx, key)
  c_vec <- numeric(n)
  for (off in .ROOK) {
    nkey <- paste(ix + off[1], iy + off[2], label, sep = ",")
    j <- lookup[nkey]                 # NA where no such neighbor exists
    hit <- !is.na(j)
    c_vec[hit] <- c_vec[hit] + y[j[hit]]
  }
  c_vec
}

# Sparse adjacency A (N x N, symmetric 0/1) for the Gibbs predictor, where the
# neighbor term depends on the UPDATED y_new and so genuinely needs the graph.
# Reused by later milestones (CSR spatial operator).
build_adjacency <- function(coords, label = NULL) {
  n  <- nrow(coords)
  ix <- as.integer(round(coords[, 1]))
  iy <- as.integer(round(coords[, 2]))
  if (is.null(label)) label <- rep(1L, n)
  key <- paste(ix, iy, label, sep = ",")
  lookup <- setNames(seq_len(n), key)
  from <- integer(0); to <- integer(0)
  for (off in .ROOK) {
    nkey <- paste(ix + off[1], iy + off[2], label, sep = ",")
    j <- lookup[nkey]
    hit <- which(!is.na(j))
    from <- c(from, hit); to <- c(to, as.integer(j[hit]))
  }
  sparseMatrix(i = from, j = to, x = 1, dims = c(n, n))
}
