####################################
# Modification of change_FF function
####################################

# samp: (p x n) data matrix
# M: number of simulated samples

change_FF_Aue <- function(samp, M = 1000)
{
	N = ncol(samp)
	D = nrow(samp)
	Sn = (1:N)
	Sn[1] = 0
	for(j in (2:N))
	{
    	Sn[j] = sum((rowSums(samp[, 1:j]) - (j/N) * rowSums(samp[,1:N]))^2)/N
	}
	k.star = min(which(Sn == max(Sn)))
	Tn = max(Sn)
	LRC_est = long_run_covariance_estimation(samp)
	lambda = eigen(LRC_est$C_0_est)$values
	asymp <- function(N)
	{
        BridgeLam = matrix(0, D, N)
        for(j in (1:D))
        {
            BridgeLam[j, ] = lambda[j] * (BBridge(0, 0, 0, 1, N - 1)^2)
        }
        return(max(colSums(BridgeLam)))
    }
	Values = sapply(1:M, function(k) asymp(N))
	z = Tn <= Values
	p = length(z[z == TRUE])/length(z)
	return(list(pvalue = p, change = k.star))
}

