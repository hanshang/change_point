##################
# load R packages
##################

require(fds)
require(ftsa)
require(strucchange)
require(ecp)

###############
# rainbow plot
###############

fert_ages = 15:49
fert_years = 1921:2021

# read data

Fert_AUS = read.csv("Fert_AUS.csv")
AUS_fert = as.matrix(cbind(Australiafertility$y, Fert_AUS[,43:48]))
colnames(AUS_fert) = fert_years

# check stationarity

T_stationary(AUS_fert)

# save figure

savepdf("aus_fert", width = 12, height = 10, toplines = 0.8)
plot(fts(fert_ages, AUS_fert), xlab = "Age", ylab = "Fertility rate")
dev.off()

############################
# regression-based approach
############################

# data_set: (p x n) data matrix
# ini_point: initial data points

AUS_fert_fun <- function(data_set, ini_point = 3)
{
    fert_years = colnames(data_set)
    fert_err_female = matrix(NA, length(fert_ages), length(fert_years) - ini_point)
    for(iw in ini_point:(length(fert_years) - 1))
    {
        # apply NFDF model
    
        NFDF_fert_female = NFDF_model(data = data_set[,1:iw])
    
        # compute errors
    
        fert_err_female[,iw - (ini_point - 1)] = data_set[,(iw + 1)] - NFDF_fert_female$forecast_curve
        print(iw + 1); rm(iw); rm(NFDF_fert_female)
    }
    colnames(fert_err_female) = (ini_point + 1):length(fert_years)
    ISFE_fert_female = ts(colSums(fert_err_female^2), start = (ini_point + 1), end = length(fert_years))
    
    # classical method
    strucchange_detect_fert_female = fert_years[(4:length(fert_years))[max(which((breakfactor(breakpoints(ISFE_fert_female ~ 1), breaks = 1)) == "segment1"))]] # 1992
    
    # ecp
    ecp_detect_fert_female = fert_years[(4:length(fert_years))[e.divisive(X = matrix(ISFE_fert_female, ncol = 1), k = 1)$estimates[2]]] # 1991
    return(c(strucchange_detect_fert_female, ecp_detect_fert_female))
}    
    
AUS_fert_fun(data_set = AUS_fert[,1:91]) # 1992 1982

# fully functional detection method

as.numeric(colnames(AUS_fert))[change_FF_Aue(samp = AUS_fert)$change] # 1976

###############################
# Booth-Maindonald-Smith (BMS)
###############################

AUS_fert_train = AUS_fert[,1:91]
AUS_fert_test = AUS_fert[,92:101]
AUS_fert_decomp = ftsm(fts(fert_ages, AUS_fert_train), order = 1)
kt = AUS_fert_decomp$coeff[,2]
m = length((fert_years)[1:91])
x <- 1:m
bp <- strucchange::breakpoints(kt ~ x)$breakpoints
minperiod = 20
bp <- bp[bp <= (m - minperiod)]
bestbreak <- max(bp)
year[(bestbreak + 1):m] # 1974 (fitting period from 1975 to 2011)
           
####################
# forecast accuracy
####################

as.numeric(colnames(AUS_fert_train))[change_FF_Aue(samp = AUS_fert[,1:91])$change] # 1974

AUS_fert_forecast_full = AUS_fert_forecast_changeFF = 
AUS_fert_forecast_strucchange = AUS_fert_forecast_bms = matrix(NA, length(15:49), 10)
for(iw in 1:10)
{
    AUS_fert_forecast_full[,iw] = forecast(ftsm(fts(15:49, AUS_fert[,1:(90+iw)]), order = 1), h = 1)$mean$y
    AUS_fert_forecast_changeFF[,iw] = forecast(ftsm(fts(15:49, AUS_fert[,54:(90+iw)]), order = 1), h = 1)$mean$y
    AUS_fert_forecast_strucchange[,iw] = forecast(ftsm(fts(15:49, AUS_fert[,72:(90+iw)]), order = 1), h = 1)$mean$y
    AUS_fert_forecast_bms[,iw] = forecast(ftsm(fts(15:49, AUS_fert[,53:(90+iw)]), order = 1), h = 1)$mean$y
    rm(iw)
}

# save figures

savepdf("AUS_fert_forecast_1", width = 12, height = 10, toplines = 0.8)
plot(fts(fert_ages, AUS_fert_test), ylim = c(0, 140), xlab = "", ylab = "Fertility rates from 2012 to 2021")
dev.off()

savepdf("AUS_fert_forecast_2", width = 12, height = 10, toplines = 0.8)
plot(fts(fert_ages,AUS_fert_forecast_full), ylim = c(0, 140), xlab = "", ylab = "Fitting period from 1921")
dev.off()

savepdf("AUS_fert_forecast_3", width = 12, height = 10, toplines = 0.8)
plot(fts(fert_ages,AUS_fert_forecast_changeFF), ylim = c(0, 140), xlab = "Age", ylab = "Fitting period from 1974")
dev.off()

savepdf("AUS_fert_forecast_4", width = 12, height = 10, toplines = 0.8)
plot(fts(fert_ages,AUS_fert_forecast_strucchange), ylim = c(0, 140), xlab = "Age", ylab = "Fitting period from 1992")
dev.off()

# MAPE forecast error

round(ftsa:::mape(forecast = AUS_fert_forecast_full, true = AUS_fert_test), 4)        # 22.6887
round(ftsa:::mape(forecast = AUS_fert_forecast_changeFF, true = AUS_fert_test), 4)    # 17.4193
round(ftsa:::mape(forecast = AUS_fert_forecast_strucchange, true = AUS_fert_test), 4) # 15.3746
round(ftsa:::mape(forecast = AUS_fert_forecast_bms, true = AUS_fert_test), 4)         # 17.4712

