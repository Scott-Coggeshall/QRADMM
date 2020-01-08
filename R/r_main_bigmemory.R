r_main_bigmemory <- function(){






}


worker_bigmemory <- function(dat_desc, inverse_desc, constraint_desc,
 beta_desc, eta_desc, i){

  dat_indices <- if(i != M) ((i - 1)*M + 1):(i*M) else ((i - 1)*M + 1):n
  inverse_indices <- ((i - 1)*p + 1):(i*p)
  dat <- bigmemory::attach.big.matrix(bm_desc)
  beta_mat <- bigmemory::attach.big.matrix(beta_desc)
  eta_mat <- bigmemory::attach.big.matrix(eta_desc)
  inverse <- bigmemory::attach.big.matrix(inverse_desc)

  constraints <- bigmemory::attach.big.matrix(constraint_desc)

  xbeta <- alpha*dat[dat_indices,2:ncol(dat) ]%*%beta_mat[, i] +
            (1 - alpha)*(dat[dat_indices, 1] - constraints[dat_indices, 1])

  constraints[dat_indices, 1] <- shrink(constraints[dat_indices, 2]/rhon +

  dat[dat_indices,1] - xbeta - .5*(2*tau - 1)/(n*rhon),
    .5*rep(1, length(dat_indices))/(n*rhon))

  beta_mat[, i] <- inverse[inverse_indices,]%*%(t(dat[dat_indices, -1])%*%dat[dat_indices, 1] - constraints[dat_indices, 1] +
    constraints[dat_indices, 2]/rhon) - eta_mat[, i]/rhon + beta_global_i)

  constraints[dat_indices, 2] <- constraints[dat_indices, 2] + rhon*(dat[dat_indices, 1] - xbeta - constraints[dat_indices, 1])

  eta_mat[, i] <- eta_mat[, i] + rhon*(beta_mat[, i] - beta_global_i)

  NULL

}
