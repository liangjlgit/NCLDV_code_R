transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- df[,1]
  rownames(t_df) <- colnames(df)
  t_df <- t_df[-1,]
  return(t_df)
}