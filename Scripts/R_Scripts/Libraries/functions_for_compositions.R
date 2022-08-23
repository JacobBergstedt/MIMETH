

balance_preds <- function(df, sbp) {
  

  
  # sbp <- read_xlsx("./Data/Cell_subset_partitions.xlsx")
  sbp <- sbp[-1]
  sbp <- as.data.frame(sbp)
  balance_basis <- balances(df, sbp)
  out <- as.data.frame(balance_basis$balances)
  names(out) <- names(sbp)
  out
  
}
