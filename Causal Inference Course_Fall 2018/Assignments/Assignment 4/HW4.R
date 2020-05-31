
## Q1: load the data and choose confounders (Step 1)
library(dplyr)
library(arm)
library(MatchIt)

# load the data
load("hw4.Rdata")
str(hw4)

# select children with birth weight less than 3000 grams
hw4_bw <- hw4[hw4$bw < 3000,]
summary(hw4_bw)

# explore the data
hw4_bw %>% group_by(treat) %>% 
          summarise(n_children = n(),
                    ppvtr36_avg = mean(ppvtr.36),
                    ppvtr36_sd = sd(ppvtr.36)/sqrt(n_children))

t.test(ppvtr.36 ~ treat, data = hw4_bw)

table(sum(hw4_bw$hispanic, hw4_bw$black, hw4_bw$white)) # the three race variables

# create a new mom education binary variable: some college or more= 1, highschool or less = 0
hw4_bw$momed_recode <- ifelse(hw4_bw$momed >2, 1, 0)
table(hw4_bw$momed_recode, hw4_bw$momed, useNA = "ifany")

# select confounder variables and create a new data set for use
confounders <- c("momage","momed_recode","prenatal","cig","booze","first","bw","preterm","black","hispanic")
new <- hw4_bw[,c("ppvtr.36", "treat", confounders)]

 # check the covariates
new %>% group_by(treat) %>% dplyr::select(one_of(confounders)) %>%
       summarise_all(funs(mean(., na.rm = T)))

#lapply(confounders, function(x) t.test(hw4_bw[, x] ~ hw4_bw[, 'treat']))

#-------------------------------------------------------------------------------
## Q2: Estimate the propensity score (Step 2)
 # use logistic regression to estimate the propensity score
log.fit <- glm (treat ~ ., data = new[, -1], family = binomial(link = "logit"))
summary(log.fit)
# save the predicted propensity scores
pscore_log <- data.frame (pscores = predict(log.fit, type = "response"),
                         treat = log.fit$model$treat)

#----------------------------------------------------------------
## Q3: Restructure your data through matching [Step 3]
 # one-to-one nearest neighbor matching with replacement
match_log <- matching(z = new$treat, score = pscore_log$pscores, replace = T)

# create "weight" variable
wts_log <- data.frame(treat = new$treat, weight = NA)

wts_log[wts_log$treat == 0, ]$weight <- match_log$cnts
wts_log[wts_log$treat == 1, ]$weight <- 1

table(wts_log$weight, wts_log$treat, useNA = "ifany")

#-----------------------------------------------------------
### Q4: Question 4: Check overlap and balance. (Step 4)
## Q4(a)examining the overlap by propensity scores
hist(pscore_log[pscore_log$treat == 0, ]$pscores, border = "blue", main = "P-scores(logistic) by group",
     xlab = "pscores (blue - control, red - treated)")
hist(pscore_log[pscore_log$treat == 1, ]$pscores, border = "red", add =T)

# examining the overlap by mother's age (momage)
hist(new[new$treat == 0,]$momage, border = "blue", main = "Mother's age by group", 
     xlab = "monther's age (blue - control, red - treated)", xlim = c(10, 45))
hist(new[new$treat == 1,]$momage, border = "red", add = T)

# examining the overlap by preterm
hist(new[new$treat == 0,]$preterm, border = "blue", main = "Preterm by group",
     xlab = "birth weight (blue - control, red - treated)")
hist(new[new$treat == 1,]$preterm, border = "red", add = T)

## Q4(c) Examining balance
# my balance checking function
check_balance <- function(df, covs, wts_df){
  
      treated <- df[df$treat == 1, ]
      control <- df[df$treat == 0, ]
      
      treated_wts<- wts_df[wts_df$treat == 1, "weight"]
      control_wts<- wts_df[wts_df$treat == 0, "weight"]
      
      bi_var <- ifelse(sapply(covs, function(x) length(unique(df[, x]))) == 2, "Y", "N") # identify binary variable
      weighted_sd <- function(x, w) sqrt(sum(w*(x - weighted.mean(x, w))^2)/sum(w)) # weighted sd function
      
      mn1 <- sapply(covs, function(x) mean(treated[, x]))
      mn0 <- sapply(covs, function(x) mean(control[, x]))
       
      mn1.m <- sapply(covs, function(x) weighted.mean(treated[, x], treated_wts))
      mn0.m <- sapply(covs, function(x) weighted.mean(control[, x], control_wts))
  
      diff <- ifelse(bi_var == "Y", mn1-mn0, (mn1-mn0)/sapply(covs, function(x) sd(treated[, x])))
      diff.m <- ifelse(bi_var == "Y", mn1.m-mn0.m, (mn1.m-mn0.m)/sapply(covs, function(x) sd(treated[, x])))
       
      ratio <- ifelse(bi_var == "N", sapply(covs, function(x) sd(control[, x])) / sapply(covs, function(x) sd(treated[, x])), NA)
      ratio.m <- ifelse(bi_var == "N", sapply(covs, function(x) weighted_sd(control[, x], control_wts)) / sapply(covs, function(x) weighted_sd(treated[, x], treated_wts)), NA)
       
      data.frame(mn1, mn0, mn1.m,mn0.m, diff, diff.m, ratio, ratio.m)
}

(balance_log <- round(check_balance(new, confounders, wts_log),2))

## Q4(e) Unit test
covs_unit <- c("bw", "b.marr")
unit_df <- hw4_bw[,c("ppvtr.36", "treat", covs_unit)]
unit_pscore <- predict(glm (treat ~ bw + b.marr, data = unit_df, family = binomial(link = "logit"))) 
unit_match <- matching(z = unit_df$treat, score = unit_pscore, replace = T)

wts_unit <- data.frame(treat = unit_df$treat, weight = NA)
wts_unit[wts_unit$treat == 0,"weight"] <- unit_match$cnts
wts_unit[wts_unit$treat == 1, "weight"] <- 1

round(check_balance(unit_df, covs_unit, wts_unit),3)

#----------------------------------------
##Q5-01 - probit regression + k1 no replacement matching
probit.fit <- glm(treat ~ ., data = new[,-1], family = binomial(link = "probit"))
summary(probit.fit)

pscore_probit <- data.frame(pscores = predict(probit.fit, type = "response"), treat = probit.fit$model$treat)
match_probit <- matching(z=new$treat, score = pscore_probit$pscores, replace = F)

wts_probit <- data.frame(treat = new$treat, matched = match_probit$matched,weight = NA)
wts_probit$weight <- ifelse(wts_probit$matched != 0, 1, 0)
table(treat = wts_probit$treat, weight = wts_probit$weight, useNA = "ifany")

hist(pscore_probit[pscore_probit$treat == 0, ]$pscores, border = "blue", main = "P-scores(probit) by group",
     xlab = "pscores (blue - control, red - treated)")
hist(pscore_probit[pscore_probit$treat == 1, ]$pscores, border = "red", add =T)

(balance_probit <- round(check_balance(new, confounders, wts_probit),2))

#------------

##Q5-02 - Optimal matching using logit distance, no drop with less predictors
match_optimal <- matchit(treat ~ momage + cig + booze + first + black + hispanic,
                         data = new, method = "optimal", distance = "logit")

pscore_opt <- data.frame(pscores = match_optimal$dist, treat = match_optimal$treat)
wts_opt <- data.frame(treat = match_optimal$treat, weight = match_optimal$weights)
table(wts_opt$weight, wts_opt$treat, useNA = "ifany")

hist(pscore_opt[pscore_opt$treat == 0, ]$pscores, border = "blue", main = "P-scores(optimal) by group",
     xlab = "pscores (blue - control, red - treated)")
hist(pscore_opt[pscore_opt$treat == 1, ]$pscores, border = "red", add =T)

(balance_opt <- round(check_balance(new, confounders, wts_opt),2))

#---------------------------------------------
# generate new data set by removing observations with momage out of the control group
c_age_min <- min(new[new$treat == 0,]$momage)
c_age_max <- max(new[new$treat == 0,]$momage)

new_age <- new %>% filter(momage >= c_age_min, momage <= c_age_max)
tapply(new_age$momage, new_age$treat, summary)
table(new_age$treat)

###Q5-03 - GAM model using new_age dataset + K1 replacement matching 
library(mgcv)
gam.fit<- gam(treat ~ momage + cig + booze + first + bw:preterm + black + hispanic, data = new_age)
summary(gam.fit)

pscore_gam <- data.frame (pscores = predict(gam.fit, type = "response"), treat = gam.fit$model$treat)
match_gam <- matching(z = new_age$treat, score = pscore_gam$pscores, replace = T)

wts_gam <- data.frame(treat = new_age$treat, weight = NA)
wts_gam[wts_gam$treat == 0, ]$weight <- match_gam$cnts
wts_gam[wts_gam$treat == 1, ]$weight <- 1

table(wts_gam$weight, wts_gam$treat, useNA = "ifany")

hist(pscore_gam[pscore_gam$treat == 0, ]$pscores, border = "blue", main = "P-scores(GAM) by group",
     xlab = "pscores (blue - control, red - treated)")
hist(pscore_gam[pscore_gam$treat == 1, ]$pscores, border = "red", add =T)
    
(balance_gam <- round(check_balance(new_age, confounders, wts_gam), 2))  

#----------------------------------------
## Q6: IPTW 
# logistic propensity score model
log.fit2 <- glm (treat ~ momage:hispanic + cig + booze + first + bw + preterm + black,
                 data = new, family = binomial(link = "logit"))
summary(log.fit2)
pscore_log2 <- data.frame(pscores = predict(log.fit2, type = "response"), treat = log.fit2$model$treat)

# creating iptw weights for ATT
wts_iptw <- pscore_log2
wts_iptw$weight_raw <- ifelse(wts_iptw$treat == 0, wts_iptw$pscores / (1 - wts_iptw$pscores), 1)
wts_iptw[wts_iptw$weight_raw > 50,]$weight_raw <- 50 # adjust extreme weights

# normalize the weights
t_sum <- sum(wts_iptw[wts_iptw$treat == 1,]$weight_raw)
c_sum <- sum(wts_iptw[wts_iptw$treat == 0,]$weight_raw)
wts_iptw$weight <- ifelse(wts_iptw$treat == 0, wts_iptw$weight_raw*t_sum/c_sum, wts_iptw$weight_raw)

tapply(wts_iptw$weight_raw, wts_iptw$treat, summary)
tapply(wts_iptw$weight_raw, wts_iptw$treat, sum)

tapply(wts_iptw$weight, wts_iptw$treat, summary)
tapply(wts_iptw$weight, wts_iptw$treat, sum)

(balance_iptw <- round(check_balance(new, confounders, wts_iptw),2))


##Q7: Comparative balance table
balance_table <- cbind(method1 = balance_log[,c(6,8)],
                       method2 = balance_probit[,c(6,8)],
                       method3 = balance_opt[,c(6,8)],
                       method4 = balance_gam[,c(6,8)],
                       method5 = balance_iptw[,c(6,8)])


##Q8: Estimate the ATT
m1.data <- cbind(new, pscores = pscore_log$pscores, weight = wts_log$weight)
m2.data <- cbind(new, pscores = pscore_probit$pscores, weight = wts_probit$weight)
m3.data <- cbind(new, pscores = pscore_opt$pscores, weight = wts_opt$weight)
m4.data <- cbind(new_age, pscores = pscore_gam$pscores, weight = wts_gam$weight)
m5.data <- cbind(new, pscores = pscore_log2$pscores, weight = wts_iptw$weight)

summary(lm(ppvtr.36 ~ factor(treat) + weight, data = m1.data))
summary(lm(ppvtr.36 ~ factor(treat) + weight, data = m2.data))
summary(lm(ppvtr.36 ~ factor(treat) + weight, data = m3.data))
summary(lm(ppvtr.36 ~ factor(treat) + weight, data = m4.data))
summary(lm(ppvtr.36 ~ factor(treat) + weight, data = m5.data))

##Q11: Linear regression
lm.fit <- lm(ppvtr.36 ~ ., data = new)
summary(lm.fit)
