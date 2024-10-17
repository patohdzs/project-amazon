library(tidyverse)
library(stargazer)
library(ggplot2)
library(broom)
library(fixest)

load("data/calibration/carbon_accumulation.Rdata")

# Fit regression on age dummies
model_1 <- lm(ratio ~ 0 + factor(age), data = df)

# Fit regression interacting with last pasture quality
model_2 <- df %>%
  filter(last_pq != 0) %>%
  lm(ratio ~ 0 + factor(age):factor(last_pq), data = .)

stargazer(model_1, model_2, out = "plots/gamma_alpha/table_1.tex")
stargazer(model_1, model_2, out = "plots/gamma_alpha/table_1.txt")

# fit regression on theoretical law of motion
df <- df %>%
  mutate(Tp = 1 - exp(-0.045 * age))

model_2 <- lm(ratio ~ 0 + Tp, data = df)
model_3 <- lm(ratio ~ Tp, data = df)
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.tex")
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.txt")

# Fit regression interacting with mean percipitation
model_3 <- df %>%
  lm(ratio ~ 0 + factor(age) + factor(age):percip, data = .)

# Fit regression interacting with last pasture quality
model_4 <- df %>%
  lm(ratio ~ 0 + percip + factor(age), data = .)

stargazer(model_3, model_4, out = "plots/gamma_alpha/table_3.tex")
stargazer(model_3, model_4, out = "plots/gamma_alpha/table_3.txt")

# Fit regression interacting with GHI
model_5 <- df %>%
  lm(ratio ~ 0 + factor(age) + factor(age):ghi, data = .)

# Fit regression with GHI levels
model_6 <- df %>%
  lm(ratio ~ 0 + ghi + factor(age), data = .)

# Fit regression interacting with both GHI and mean percipitation
model_7 <- df %>%
  lm(ratio ~ 0 + percip + ghi + factor(age), data = .)

stargazer(model_5, model_6, model_7, out = "plots/gamma_alpha/table_4.tex")
stargazer(model_5, model_6, model_7, out = "plots/gamma_alpha/table_4.txt")

fs <- df %>%
  lm(gamma ~ percip + ghi, data = .)

df <- df |> mutate(res = fs$residuals)
