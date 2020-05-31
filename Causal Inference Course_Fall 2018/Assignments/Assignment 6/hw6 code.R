#Q1
set.seed(123)
income <- c(
  runif(1000*0.28, 0, 2), 
  runif(1000*0.07, 18, 20),
  rnorm(1000*0.30, mean=15, sd=5),
  rnorm(1000*0.30, mean=28, sd=8),
  rnorm(1000*0.05, mean=40, sd=10)
)

negs <- length(income[income<0])
income[income<0] <- runif(negs,0,60) #negative income into 0
income[income>60] <- 60 # maximal is 60

summary(income)
hist(income, xlim = c(0, 60), ylim = c(0, 0.15),breaks = 30, freq = F, main = "Simulated income (unit: $1000)")

# Q2
eligible <- ifelse(income < 20, 1, 0)
table(eligible, useNA = "ifany")

# Q3a: linear models for both y1 and y0
fullA <- data.frame(income, y0 = NA, y1 = NA)
set.seed(31)

fullA$y0 <-  5.5 + 0.25*income + rnorm(1000, mean = 0, sd = 2)
summary(fullA$y0)

fullA$y1 <- fullA$y0 + rnorm(1000, mean = 4, sd = 2)
summary(fullA$y1)

mean(fullA$y1 - fullA$y0)

library(dplyr)
obsA <- fullA %>% cbind(eligible) %>%
        mutate(y = ifelse(eligible == 1, y1, y0))%>% 
        select(income, eligible, y)

# Q3b
fullB <- data.frame(income, y0 = NA, y1 = NA)
income_new <- income - 20

set.seed(32)
fullB$y0 <- 14 + 0.4*income_new + rnorm(1000, 0, 2)
summary(fullB$y0)

for (i in 1:1000){
if(abs(income_new[i]) < 2) {
  fullB$y1[i] = fullB$y0[i] + rnorm(1, 4, 2)
} else {
  fullB$y1[i] = 1/125*income_new[i]^2 + rnorm(1, 14, 2)
}
}
summary(fullB$y1)

fullB_sub <- subset(fullB, abs(fullB$income-20) < 2)
mean(fullB_sub$y1-fullB_sub$y0)

obsB <- fullB %>% cbind(eligible) %>%
        mutate(y = ifelse(eligible == 1, y1, y0))%>%
        select(income, eligible, y)

# Q4
library(ggplot2)
ggplot(data = obsA, aes(income, y, color = as.factor(eligible))) + geom_point() +
  scale_color_manual(values = c("blue","red")) + ggtitle("Scatterplot of obsA")

ggplot(data = obsB, aes(income, y, color = as.factor(eligible))) + geom_point() +
  scale_color_manual(values = c("blue","red")) + ggtitle("Scatterplot of obsB")

#Q5a
(for.estimate <- data.frame(income = c(20, 20), eligible = c(1, 0)))

fitA.1 <- lm(y ~ income + eligible, data = obsA)
pred1 <- predict(fitA.1, for.estimate, type = "response")
pred1[[1]]-pred1[[2]]

fitA.2 <- lm(y ~ eligible * income, data = obsA)
pred2 <- predict(fitA.2, for.estimate, type = "response")
pred2[[1]]-pred2[[2]]

fitA.3 <- lm(y ~ I(income^2) * eligible * income, data = obsA)
pred3 <- predict(fitA.3, for.estimate, type = "response")
pred3[[1]]-pred3[[2]]

#Q5b
fitB.1 <- lm(y ~ income + eligible, data = obsB)
pred4 <- predict(fitB.1, for.estimate, type = "response")
pred4[[1]]-pred4[[2]]

fitB.2 <- lm(y ~ eligible * income, data = obsB)
pred5 <- predict(fitB.2, for.estimate, type = "response")
pred5[[1]]-pred5[[2]]

fitB.3 <- lm(y ~ I(income^2) * eligible * income, data = obsB)
pred6 <- predict(fitB.3, for.estimate, type = "response")
pred6[[1]]-pred6[[2]]

#Q6a
selected <- (income <= 22 & income >= 18)

rd.fitA.1 <- lm(y ~ income + eligible, data = obsA, subset = selected)
pred7 <- predict(rd.fitA.1, for.estimate, type = "response")
pred7[[1]]-pred7[[2]]

rd.fitA.2 <- lm(y ~ eligible * income, data = obsA, subset = selected)
pred8 <- predict(rd.fitA.2, for.estimate, type = "response")
pred8[[1]]-pred8[[2]]

rd.fitA.3 <- lm(y ~ I(income^2) * eligible * income, data = obsA, subset = selected )
pred9 <- predict(rd.fitA.3, for.estimate, type = "response")
pred9[[1]]-pred9[[2]]

#Q6b
rd.fitB.1 <- lm(y ~ income + eligible, data = obsB, subset = selected)
pred10 <- predict(rd.fitB.1, for.estimate, type = "response")
pred10[[1]]-pred10[[2]]

rd.fitB.2 <- lm(y ~ eligible * income, data = obsB, subset = selected)
pred11 <- predict(rd.fitB.2, for.estimate, type = "response")
pred11[[1]]-pred11[[2]]

rd.fitB.3 <- lm(y ~ I(income^2) * eligible * income, data = obsB, subset = selected )
pred12 <- predict(rd.fitB.3, for.estimate, type = "response")
pred12[[1]]-pred12[[2]]

#Q7
estimate_df <- as.data.frame(rbind(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9, pred10, pred11, pred12))
names(estimate_df) <- c("eligible", "not_eligible")
estimate_df$estimation <- estimate_df$eligible - estimate_df$not_eligible
print(estimate_df)

#Q8e
hist(income, breaks = 33)

# Challenge problem
library(rdrobust)

rdplot(obsA$y, obsA$income, c = 20, binselect = "esmv", x.lim = c(0, 60), y.lim = c(0, 30), title = "RDPlot of obsA")
rdplot(obsB$y, obsB$income, c = 20, binselect = "esmv", x.lim = c(0, 60), y.lim = c(0, 30), title = "RDPlot of obsB")

# Way 1:  MSE-optimal method (same bandwith on both sides of the cutoff)
rd.bw1 <- rdbwselect(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 1, bwselect = "mserd")
summary(rd.bw1)
# RD estimation
summary(rdrobust(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 1, bwselect = "mserd"))

# Way 2: MSE(Mean Square Error)-optimal method - different bandwidth on each side of the cutoff
rd.bw2 <- rdbwselect(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 1, bwselect = "msetwo")
summary(rd.bw2)

summary(rdrobust(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 2, bwselect = "msetwo"))

# Way 3: CER(Coverage Error Rate)-optimal-different bandwidth on each side of the cutoff
rd.bw3 <- rdbwselect(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 1, bwselect = "cercomb2")
summary(rd.bw3)

summary(rdrobust(obsB$y, obsB$income, kernel = "triangular", c = 20, p = 2, bwselect = "cercomb2"))
