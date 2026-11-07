# function to rescale colours 
scale_fill_ratio <- function(range, transform = "identity",
                             palette = "RdBu", n = 1024, ...) {
  tr  <- scales::as.transform(transform)   # scales < 1.3: as.trans()
  tf  <- tr$transform
  inv <- tr$inverse
  
  trans_range  <- tf(range)  # transformed range
  half_range   <- max(abs(trans_range))  # symmetric half-width about 0
  col_vals     <- seq(trans_range[1], trans_range[2], length.out = n) # even in transformed space
  
  pal  <- scales::gradient_n_pal(hcl.colors(256, palette, rev = TRUE))
  cols <- pal((col_vals + half_range) / (2 * half_range) )  # -> [0, 1], locked at 0 = midpoint
  
  scale_fill_gradientn(colours = cols, limits = range,
                       transform = transform, ...)
}