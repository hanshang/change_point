####################
# independence test
####################

# data: (p x n) data matrix
# K: number of retained functional principal components
# H: tuning parameter

GK_test <- function(data, K, H = 10)
{
    p = nrow(data)
    n = ncol(data)
    
    cen_data = t(scale(t(data), scale = FALSE))
    
    pc_decomp = ftsm(fts(1:p, cen_data), order = K, mean = FALSE)
    score = pc_decomp$coeff
    
    C0 = crossprod(score)/n
    
    c_h = array(NA, dim = c(K, K, H))
    for(h in 1:H)
    {
      for(k in 1:K)
      {
        for(l in 1:K)
        {
          score_uni = 0
          for(t in 1:(n - h))
          {
            score_uni = score_uni + (score[t,k] * score[t+h,l])
          }
          c_h[k,l,h] = score_uni/n
        }
      }
    }
    
    r_f_h = r_b_h = array(NA, dim = c(K, K, H))
    for(h in 1:H)
    {
      r_f_h[,,h] = solve(C0) %*% c_h[,,h]
      r_b_h[,,h] = c_h[,,h] %*% solve(C0)
    }
    
    summand = vector("numeric", H)
    for(h in 1:H)
    {
      summand[h] = sum(r_f_h[,,h] * r_b_h[,,h])
    }
    Q_n = n * sum(summand)
    
    p_value = 1 - pchisq(Q_n, df = K^2 * H)
    return(p_value)    
}
