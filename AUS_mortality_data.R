##########################
# load existing functions
##########################

source("change_FF_Aue.R")
source("long_run_covariance_h.R")
source("load_packages.R")

##################################
# Australian mortality & exposure
##################################

AUS_rate = readHMDweb(CNTRY = "AUS", item = "Mx_1x1", username = "hanlin.shang@gmail.com",
                      password = "Toyota1985@", fixup = FALSE)
AUS_expo = readHMDweb(CNTRY = "AUS", item = "Exposures_1x1", username = "hanlin.shang@gmail.com",
                      password = "Toyota1985@", fixup = FALSE)

year = 1921:2020
age = 0:110
n_year = length(year)
n_age = length(age)
AUS_rate_female = matrix(AUS_rate$Female, n_age, n_year)
AUS_rate_male   = matrix(AUS_rate$Male, n_age, n_year)

AUS_expo_female = matrix(AUS_expo$Female, n_age, n_year)
AUS_expo_male   = matrix(AUS_expo$Male, n_age, n_year)

##################
# demogdata class
##################

AUS_demo_female = demogdata(AUS_rate_female, AUS_expo_female, ages = 0:110, years = 1921:2020,
                            "mortality", "AUS", "Female")
AUS_demo_male = demogdata(AUS_rate_male, AUS_expo_male, ages = 0:110, years = 1921:2020,
                          "mortality", "AUS", "Male")

# truncating from 0 to 100+

AUS_demo_female_trunc = extract.ages(data = AUS_demo_female, ages = 0:100)
AUS_demo_male_trunc   = extract.ages(data = AUS_demo_male, ages = 0:100)

# applying the BMS method

BMS_change_point_female_mort = bms(AUS_demo_female_trunc) # 1981 (fitting period is 1982-2020)
BMS_change_point_male_mort   = bms(AUS_demo_male_trunc)   # 1973 (ftting period is 1974-2020)

#####################
# graphical displays
#####################

savepdf("AUS_female_rainbow", width = 12, height = 10, toplines = 0.8)
plot(fts(0:100, log(AUS_demo_female_trunc$rate$Female, base = 10)), xlab = "Age", ylab = "Log death rate",
     main = "AUS: female death rates (1921-2020)")
dev.off()

savepdf("AUS_male_rainbow", width = 12, height = 10, toplines = 0.8)
plot(fts(0:100, log(AUS_demo_male_trunc$rate$Male, base = 10)), xlab = "Age", ylab = "Log death rate",
     main = "AUS: male death rates (1921-2020)")
dev.off()

###########
# fd class
###########

AUS_demo_female_fd = log(AUS_demo_female_trunc$rate$Female, base = 10)
AUS_demo_male_fd   = log(AUS_demo_male_trunc$rate$Male,     base = 10)

# detect change point via change_FF function

change_FF_AUS_female = change_FF_Aue(samp = AUS_demo_female_fd)
change_FF_AUS_male   = change_FF_Aue(samp = AUS_demo_male_fd)

c(change_FF_AUS_female$pvalue, year[change_FF_AUS_female$change]) # 0.005 (subject to randomness) 1972
c(change_FF_AUS_male$pvalue, year[change_FF_AUS_male$change])     # 0.002 (subject to randomness) 1977
