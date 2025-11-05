# A statistical project in R investigating the relationship between water hardness and mortality in towns in the US and Canada
# Ethan Milne 2023

set.seed(62) # Assigned sample for project

# Folders for repository
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

# Use already-saved samples in /data
x <- read.csv("data/x_canada.csv")[[1]]
y <- read.csv("data/y_usa.csv")[[1]]
m1 <- read.csv("data/mortality_calcium_25towns.csv")

# Part 1: Compare mean calcium (Canada vs USA) ---------------------------
# Visual checks (optional in headless runs)
# qqnorm(x); qqline(x); qqnorm(y); qqline(y)

t_welch <- t.test(x, y, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
mean_diff <- mean(x) - mean(y)

# Nonparametric sensitivity
wilx <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided", conf.int = TRUE))

# --- Part B: Mortality vs Calcium (n = 25 towns) ----------------------------
mortality <- m1[[1]]
calcium   <- m1[[2]]

# Scatter with regression line
png("results/figures/scatter_mortality_calcium.png", width = 900, height = 700, res = 120)
plot(calcium, mortality,
     xlab = "Calcium in drinking water (ppm)",
     ylab = "Annual mortality (per 100,000)",
     main = "Mortality vs Calcium (25 towns)")
abline(lm(mortality ~ calcium), lty = 2)
dev.off()

# Correlation (Pearson and Spearman)
cor_pearson  <- cor.test(mortality, calcium, method = "pearson")
cor_spearman <- cor.test(mortality, calcium, method = "spearman", exact = FALSE)

# Linear model: mortality as outcome, calcium as predictor
fit <- lm(mortality ~ calcium)
pred40 <- predict(fit, newdata = data.frame(calcium = 40),
                  interval = "prediction", level = 0.95)

# Data plot
png("results/figures/lm_diagnostics.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2,2)); plot(fit); par(mfrow = c(1,1))
dev.off()

# Text report
sink("results/summary.txt")
cat("Calcium comparison (Canada vs USA)\n")
cat("==================================\n")
cat(sprintf("Mean(Canada) = %.3f ppm\n", mean(x)))
cat(sprintf("Mean(USA)    = %.3f ppm\n", mean(y)))
cat(sprintf("Mean difference (Canada - USA) = %.3f ppm\n", mean_diff))
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
            t_welch$statistic, t_welch$parameter, t_welch$p.value))
cat(sprintf("95%% CI for mean difference: [%.3f, %.3f]\n\n",
            t_welch$conf.int[1], t_welch$conf.int[2]))
cat(sprintf("Wilcoxon rank-sum (robustness): W = %.0f, p = %.4f\n\n",
            wilx$statistic, wilx$p.value))

cat("Mortality vs Calcium\n")
cat("====================\n")
cat(sprintf("Pearson r = %.3f (p = %.4f)\n", cor_pearson$estimate, cor_pearson$p.value))
cat(sprintf("Spearman rho = %.3f (p = %.4f)\n\n",
            cor_spearman$estimate, cor_spearman$p.value))

cat("Linear model: mortality ~ calcium\n")
cat("---------------------------------\n")
print(summary(fit))
cat("\nPrediction for calcium = 40 ppm (95% prediction interval):\n")
print(pred40)
cat("\nSession info:\n")
print(sessionInfo())
sink()
