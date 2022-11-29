CoDF <- function(da, daa) {
  up <- upper.tri(da)
  data.frame(
    from = rownames(da)[col(da)[up]],
    to = rownames(da)[row(da)[up]],
    corr  =(da)[up],
    p = daa[up]
  )
}