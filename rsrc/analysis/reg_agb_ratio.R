library(ggplot2)

load("data/calibration/combined_df.Rdata")

combined_df_pq <- combined_df %>%
  filter(!is.na(agb) & !is.na(sec) & !is.na(gamma) & pq !=0) %>%
  mutate(ratio = agb / gamma)

combined_df <- combined_df %>%
  filter(!is.na(agb) & !is.na(sec) & !is.na(gamma)) %>%
  mutate(ratio = agb / gamma)



generate_dummies <- function(df, sec_var, max_sec) {
  for (i in 1:max_sec) {
    dummy_name <- paste0("dummy_", i)
    df <- df %>%
      mutate(!!dummy_name := ifelse(sec > (i - 1) & sec <= i, 1, 0))
  }
  return(df)
}

# Apply the function to generate 32 dummy variables
combined_df <- generate_dummies(combined_df, "sec", 32)


formula <- as.formula(paste("ratio ~  ", paste(paste0("dummy_", 2:32), collapse = " + ")))


model <- lm(formula , data = combined_df)
summary(model)



coefficients <- coef(model)[paste0("dummy_", 2:32)]
stderr <- summary(model)$coefficients[paste0("dummy_", 2:32), "Std. Error"]

coefficients_df <- data.frame(
  dummy = 2:32,
  coefficient = coefficients,
  lower_bound = coefficients - stderr,
  upper_bound = coefficients + stderr,
  type = "Coefficients"
)

# Generate the theoretical function values starting from t=0 to t=32
x_values <- 0:32
theoretical_values <- 1 - exp(-0.045 * x_values)

# Create a data frame for plotting the theoretical function
theoretical_df <- data.frame(
  dummy = x_values,
  coefficient = theoretical_values,
  type = "Theoretical Function"
)

# Combine the two data frames
plot_df <- bind_rows(coefficients_df, theoretical_df)

# Plot the coefficients with error bars and overlay the theoretical function
p <- ggplot(plot_df, aes(x = dummy, y = coefficient, color = type, linetype = type)) +
  geom_point(data = coefficients_df) +
  geom_line(data = coefficients_df) +
  geom_errorbar(data = coefficients_df, aes(ymin = lower_bound, ymax = upper_bound), width = 0.2) +
  geom_line(data = theoretical_df, size = 1) +
  labs(title = "Coefficients of Dummy Variables with Theoretical Function",
       x = "Dummy Variable Index (age)",
       y = "Coefficient Value",
       color = "Legend",
       linetype = "Legend") +
  scale_color_manual(values = c("Coefficients" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Coefficients" = "solid", "Theoretical Function" = "dashed")) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.15))  # Position the legend inside the plot



# Save the plot
ggsave("coefficients_plot_with_intercept.png", plot = p, width = 8, height = 6)


stop()

combined_df <- combined_df %>%
  mutate(Tp = 1-exp(-0.045*sec))


model2 <- lm(ratio ~ -1 +  Tp , data = combined_df)
summary(model2)



