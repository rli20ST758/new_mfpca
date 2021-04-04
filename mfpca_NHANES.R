## The aim is to do MFPCA on NHANES data using existing software
## After getting MFPCA result, fitting AFCM on subject average and by visit

rm(list = ls())
library(mgcv)
source("./code/mfpca.fast.R")
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

## load data
df_MIMS_1440 <- readRDS("./data/NHANES_2011_2014_MIMS_1440.rds")
df_MIMS_features <- readRDS("./data/NHANES_2011_2014_MIMS_features.rds")

## replace MIMS = -0.01 with NA. These values indicate MIMS could not be calculated for this minute.
df_MIMS_1440[df_MIMS_1440 < 0] <- NA

# ############################################################################
# ## do I2C2 on minute-level MIMS
# ############################################################################
# 
# fit_I2C2 <- I2C2(y = df_MIMS_1440, id = df_MIMS_features$SEQN, visit = df_MIMS_features$day_of_week)
# fit_I2C2$lambda ## I2C2 = 0.195, indicating the subjects/visits are more independent than just replicate

############################################################################
## do MFPCA on minute-level MIMS
############################################################################

ptm <- proc.time()
fit_mfpca <- mfpca.fast(Y = df_MIMS_1440, id = df_MIMS_features$SEQN, 
                        group = df_MIMS_features$day_of_week, silent = FALSE)
proc.time() - ptm ## 524s!

# fit_mfpca <- readRDS("./results/mfpca_NHANES.rds")


## build the fixed dataset
dat_fixed <- as.data.frame(fit_mfpca$eta + fit_mfpca$mu)
colnames(dat_fixed) <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
dat_fixed$Overall <- fit_mfpca$mu
dat_fixed$time <- 1:1440
dat_fixed_long <- pivot_longer(dat_fixed, cols = 1:8, names_to = "Day")
dat_fixed_long$Day <- factor(dat_fixed_long$Day, levels = c("Overall", "Sunday", "Monday", 
                             "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))

## plot mu and eta
pdf("./figures/nhanes_mueta.pdf", width = 7, height = 3)
ggplot(dat_fixed_long, aes(x = time, y = value, group = Day)) +
  theme_classic() +
  geom_line(aes(linetype = Day, color = Day)) +
  scale_linetype_manual(values = c("solid", "dotdash", rep("dashed",5), "dotdash")) +
  scale_color_manual(values=c("black",rainbow(7))) +
  scale_x_continuous(breaks = c(1,6,12,18,23)*60, 
                     labels = c("01:00","06:00","12:00","18:00","23:00")) +
  labs(x = "Time of day", y = "MIMS", color = "", linetype = "") +
  theme(legend.position="right", 
        axis.title = element_text(size = 10))
dev.off()

# ## plot visit-specific mean curves
# jpeg("./figures/mfpca_eta.jpeg", width = 1200, height = 700, res = 130)
# matplot(fit_mfpca$eta, type = "l", col = 1:7, xlab = "Time of Day",
#         ylab = expression(eta[j](s)), main = "Day of week specific mean curves", xaxt = "n")
# axis(1, at=(c(1,6,12,18,23)*60), labels=c("01:00","06:00","12:00","18:00","23:00"), cex.axis=1)
# nn <- ncol(fit_mfpca$eta)
# legend("bottomright", legend = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"),
#        col=seq_len(nn),cex=1,fill=seq_len(nn))
# dev.off()

## plot level 1 first 6 eigenfunctions
# jpeg("./figures/mfpca_level1.jpeg", width = 1200, height = 700, res = 130)
pdf("./figures/nhanes_level1.pdf", width = 9, height = 2.5)
# par(mfrow = c(2, 3))
par(mfrow = c(1, 3))
for(i in 1:3){
  plot(fit_mfpca$efunctions$level1[,i], type = "l", xlab = "Time of day", xaxt = "n", ylab = "", cex.lab=1.4, cex.main = 1.4,
       main = paste0("level 1, comp ", i, ", ", round(fit_mfpca$evalues$level1[i]/sum(fit_mfpca$evalues$level1)*100, 2), "%"))
  axis(1, at=(c(1,6,12,18,23)*60), labels=c("01:00","06:00","12:00","18:00","23:00"))
  abline(h = 0, lty = 2)
}
dev.off()

## plot level 2 first 6 eigenfunctions
# jpeg("./figures/mfpca_level2.jpeg", width = 1200, height = 700, res = 130)
pdf("./figures/nhanes_level2.pdf", width = 9, height = 2.5)
# par(mfrow = c(2, 3))
par(mfrow = c(1, 3))
for(i in 1:3){
  plot(fit_mfpca$efunctions$level2[,i], type = "l", xlab = "Time of day", xaxt = "n", ylab = "", cex.lab=1.4, cex.main = 1.4,
       main = paste0("level 2, comp ", i, ", ", round(fit_mfpca$evalues$level2[i]/sum(fit_mfpca$evalues$level2)*100, 2), "%"))
  axis(1, at=(c(1,6,12,18,23)*60), labels=c("01:00","06:00","12:00","18:00","23:00"))
  abline(h = 0, lty = 2)
}
dev.off()

## variance explained by within cluster variability
rhoW <- sum(fit_mfpca$evalues$level1) / (sum(fit_mfpca$evalues$level2) + sum(fit_mfpca$evalues$level1)) ## 0.31







# ############################################################################
# ## do AFCM on subject curves Z_i(s)
# ############################################################################
# 
# ## mortality and covariates dataset
# ## Among 19931 subjects, we exclude subjects with
# ## 1. age less than 50 or greater than or equal to 80 (15162 subjects);
# ## 2. have missing variable of interest (964 subjects).
# ## The final dataset contains 3805 subjects, with 150 events before Dec 31, 2015.
# df_analysis <- readRDS("./data/NHANES_2011_2014_mort_covar.rds")
# 
# ## obtain Z_i(s)
# ind.analy <- which(unique(df_MIMS_features$SEQN) %in% df_analysis$SEQN)
# Zhat <- fit_mfpca$scores$level1[ind.analy,] %*% t(fit_mfpca$efunctions$level1)
# 
# df_final <- df_analysis
# df_final$Zhat_q <- I(apply(Zhat, 2, function(y) ecdf(y)(y))) # quantiles
# 
# ## add variables for numerical integration
# df_final$lmat <- I(matrix(1/1440, ncol = 1440, nrow = nrow(df_final)))
# df_final$tmat <- I(matrix(seq(0, 1, length.out = 1440), ncol = 1440, 
#                           nrow = nrow(df_final), byrow = TRUE))
# 
# ## fit afcm
# fit_Zhat_q <- gam(event_time_years ~ gender + age_years_interview + education_adult + race + BMI_cat +
#                   overall_health_combined + diabetes + heart_attack + CHF + CHD + stroke + cancer +
#                   mobility_problem + alcohol_consumption_fac + cigarette_smoking +
#                   ti(tmat, Zhat_q, by = lmat, bs = c("cc", "cr"), k = c(5, 5), mc = c(FALSE, TRUE)),
#                 weights = mortstat, data = df_final, family = cox.ph())
# jpeg("./figures/mfpca_afcm.jpeg", width = 1000, height = 800, res = 100)
# vis.gam(fit_Zhat_q, view=c("tmat","Zhat_q"), plot.type="contour", color="topo", xaxt = "n",
#         xlab = "Time of Day", ylab = "Quantile of Z hat", main= paste0("Estimates using subject curve"))
# axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"))
# dev.off()
# 
# ############################################################################
# ## do AFCM on subject day of week curves Z_i(s) + W_{ij}(s)
# ############################################################################
# 
# fit_ZWhat_q <- list()
# 
# # dow <- 1 ## specify day of week (1 = Sunday, ...)
# for(dow in 1:7){
# df_MIMS <- df_MIMS_features[which(df_MIMS_features$day_of_week == dow), 1:5]
# 
# ## First calculate W for subjects with visit at this day of week
# df_MIMS$W <- I(fit_mfpca$scores$level2[which(df_MIMS_features$day_of_week == dow),] %*% t(fit_mfpca$efunctions$level2))
# df_final <- left_join(df_analysis, df_MIMS, by = "SEQN")
# ## remove subjects without accelerometry data of that day of week
# df_final <- df_final[which(!is.na(df_final$day_of_week)),]
# 
# ## Second calculate Z for remaining subjects in the analysis
# ind.analy <- which(unique(df_MIMS_features$SEQN) %in% df_final$SEQN)
# df_final$Z <- I(fit_mfpca$scores$level1[ind.analy,] %*% t(fit_mfpca$efunctions$level1))
# 
# ## Finally calculate Z + W and its quantile
# df_final$ZWhat <- I(df_final$Z + df_final$W)
# df_final$ZWhat_q <- I(apply(df_final$ZWhat, 2, function(y) ecdf(y)(y))) # quantiles
# 
# ## add variables for numerical integration
# df_final$lmat <- I(matrix(1/1440, ncol = 1440, nrow = nrow(df_final)))
# df_final$tmat <- I(matrix(seq(0, 1, length.out = 1440), ncol = 1440, 
#                           nrow = nrow(df_final), byrow = TRUE))
# 
# ## fit afcm
# fit_ZWhat_q[[dow]] <- gam(event_time_years ~ gender + age_years_interview + education_adult + race + BMI_cat +
#                     overall_health_combined + diabetes + heart_attack + CHF + CHD + stroke + cancer +
#                     mobility_problem + alcohol_consumption_fac + cigarette_smoking +
#                     ti(tmat, ZWhat_q, by = lmat, bs = c("cr", "cr"), k = c(5, 5), mc = c(FALSE, TRUE)),
#                     weights = mortstat, data = df_final, family = cox.ph())
# print(dow)
# }
# 
# jpeg("./figures/mfpca_afcm.jpeg", width = 2000, height = 1000, res = 150)
# par(mfrow = c(2, 4))
# vis.gam(fit_Zhat_q, view=c("tmat","Zhat_q"), plot.type="contour", color="topo", xaxt = "n",
#         xlab = "Time of Day", ylab = "Quantile of Z hat", main= paste0("Estimates using subject curve"))
# axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"))
# 
# for(dow in 1:7){
#   vis.gam(fit_ZWhat_q[[dow]], view=c("tmat","ZWhat_q"), plot.type="contour", color="topo", xaxt = "n",
#           xlab = "Time of Day", ylab = "Quantile of Z hat", main= paste0("Estimates, Day of week = ", dow))
#   axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("00:00", "06:00", "12:00", "18:00", "24:00"))
# }
# dev.off()


