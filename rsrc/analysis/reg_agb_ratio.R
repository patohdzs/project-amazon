library(tidyverse)
library(stargazer)
library(ggplot2)

load("data/calibration/carbon_accumulation.Rdata")

# Filter out empty pixels
df <- df %>%
  filter(!is.na(co2e) & !is.na(age) & !is.na(gamma) & !is.na(pq))

# Pre-processing
df <- df %>%
  mutate(age = round(age)) %>%
  mutate(ratio = co2e / gamma)

# Fit regression from approach 1 (NOTE exclude level 1)
model_1 <- lm(ratio ~ 0 + factor(age), data = df)
stargazer(model_1, out = "plots/gamma_alpha/table_1.tex")
stargazer(model_1, out = "plots/gamma_alpha/table_1.txt")

# Get coefficients
coefs <- coef(model_1)[paste0("factor(age)", 2:32)]
coefs_df <- data.frame(
  time = 2:32,
  coef = coefs,
  type = "Coefficients"
)

# Generate the theoretical values from t=0 to 32
time <- 0:32
theoretical_values <- 1 - exp(-0.045 * time)

# Create a data frame for plotting the theoretical function
theoretical_df <- data.frame(
  time = time,
  coef = theoretical_values,
  type = "Theoretical Function"
)

# Combine the two data frames
plot_df <- bind_rows(coefs_df, theoretical_df)

# Plot the coefs and overlay the theoretical function
p <- plot_df |> ggplot(aes(x = time, y = coef, color = type, linetype = type)) +
  geom_point(data = coefs_df) +
  geom_line(data = coefs_df) +
  geom_line(data = theoretical_df, size = 1) +
  labs(
    x = "Time (years)",
    y = "Coefficient",
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c("Coefficients" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Coefficients" = "solid", "Theoretical Function" = "dashed")) +
  theme_minimal()

# Save the plot
ggsave("plots/gamma_alpha/coefs_plot.png", plot = p)

# Perform approach 2 regression
df <- df %>%
  mutate(Tp = 1 - exp(-0.045 * age))

model_2 <- lm(ratio ~ 0 + Tp, data = df)
model_3 <- lm(ratio ~ Tp, data = df)
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.tex")
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.txt")
