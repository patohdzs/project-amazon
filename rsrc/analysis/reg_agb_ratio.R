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


formula <- as.formula(paste("ratio ~ -1 +", paste(paste0("dummy_", 1:32), collapse = " + ")))


model <- lm(formula , data = combined_df)
summary(model)





coefficients <- coef(model)[paste0("dummy_", 1:32)]

# Create a data frame for plotting
coefficients_df <- data.frame(
  dummy = 1:32,
  coefficient = coefficients
)

# Plot the coefficients
p <- ggplot(coefficients_df, aes(x = dummy, y = coefficient)) +
  geom_point() +
  geom_line() +
  labs(title = "Coefficients of Dummy Variables",
       x = "Dummy Variable Index (age)",
       y = "Coefficient Value") +
  theme_minimal()
ggsave("coefficients_plot.png", plot = p, width = 8, height = 6)


combined_df <- combined_df %>%
  mutate(Tp = 1-exp(-0.045*sec))


model2 <- lm(ratio ~  Tp , data = combined_df)
summary(model2)



