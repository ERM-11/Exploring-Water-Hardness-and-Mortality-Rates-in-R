# Project: Water hardness (calcium, ppm) and mortality
# Author: Ethan Milne
# Brief: Simple analyses and clear plots that directly answer the assessed questions.
# Data expected:
#   data/x_canada.csv              (vector of calcium ppm for 20 Canadian towns)
#   data/y_usa.csv                 (vector of calcium ppm for 20 US towns)
#   data/mortality_calcium_25towns.csv  (two columns: mortality, calcium, n=25 towns)

# -------------------------- Setup --------------------------------------------
k <- 62  # last two digits of student ID as given in the report
set.seed(k)

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

need <- c("data/x_canada.csv", "data/y_usa.csv", "data/mortality_calcium_25towns.csv")
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop(
    "Missing input files:\n",
    paste0(" - ", missing, collapse = "\n"),
    "\nPlace them in the paths shown and re-run."
  )
}

# -------------------------- Load data ----------------------------------------
# Part A data: two independent samples of calcium concentration (ppm)
x_can <- read.csv("data/x_canada.csv", header = TRUE, check.names = FALSE)[[1]]
y_usa <- read.csv("data/y_usa.csv",   header = TRUE, check.names = FALSE)[[1]]

# Part B data: paired (row-wise) mortality and calcium for 25 towns
m1 <- read.csv("data/mortality_calcium_25towns.csv", header = TRUE, check.names = FALSE)
if (ncol(m1) < 2) stop("Expected two columns in mortality_calcium_25towns.csv: mortality, calcium")
mortality <- m1[[1]]
calcium   <- m1[[2]]

# -------------------------- Part A: Compare means ----------------------------
# Aim: Test if mean calcium (Canada) equals mean calcium (USA).
# Hypotheses:
#   H0: mu_CAN - mu_USA = 0
#   H1: mu_CAN - mu_USA != 0
# Test: Welch two-sample t-test (allows unequal variances), alpha = 0.05.
# Assumptions: independent random samples; each sample approximately normal.
#   We check approximate normality with Shapiro-Wilk (n=20 per group).

shap_can <- tryCatch(shapiro.test(x_can), error = function(e) NULL)
shap_usa <- tryCatch(shapiro.test(y_usa), error = function(e) NULL)

t_welch <- t.test(x_can, y_usa, alternative = "two.sided", var.equal = FALSE, conf.level = 0.95)
mean_diff <- mean(x_can) - mean(y_usa)

# Visual: simple boxplot (no extras)
png("results/figures/boxplot_calcium_by_country.png", width = 900, height = 650, res = 120)
par(mar = c(5,5,4,2) + 0.1)
boxplot(list(Canada = x_can, USA = y_usa),
        ylab = "Calcium concentration (ppm)",
        main = "Calcium in Drinking Water: Canada vs USA")
grid()
dev.off()

# -------------------------- Part B: Mortality vs Calcium ---------------------
# (a) Scatter with fitted least-squares line
png("results/figures/scatter_mortality_calcium.png", width = 900, height = 650, res = 120)
par(mar = c(5,5,4,2) + 0.1)
plot(calcium, mortality,
     xlab = "Calcium in drinking water (ppm)",
     ylab = "Annual mortality (per 100,000)",
     main = "Mortality vs Calcium (25 towns)",
     pch = 19)
abline(lm(mortality ~ calcium), lty = 2, lwd = 2)
grid()
dev.off()

# (b) Test H0: rho = 0 vs H1: rho != 0 using Pearson correlation (as per brief)
cor_pearson <- cor.test(mortality, calcium, method = "pearson")

# (c) Regression of mortality on calcium, and 95% prediction interval at calcium = 40 ppm
fit <- lm(mortality ~ calcium)
pred40 <- predict(fit, newdata = data.frame(calcium = 40),
                  interval = "prediction", level = 0.95)

# Reliability note: report whether x=40 is inside observed range
calcium_range <- range(calcium, na.rm = TRUE)
inside_range <- (40 >= calcium_range[1] && 40 <= calcium_range[2])

# -------------------------- Report (plain text) ------------------------------
sink("results/summary.txt"); on.exit(sink(), add = TRUE)

cat("Statistics Project: Water hardness (calcium) and mortality\n")
cat("Author: Ethan Milne   |   k =", k, "\n\n")

cat("Q1. How to obtain random samples (concise plan)\n")
cat(" - Define sampling frame: list of all large towns in Canada and in the USA.\n")
cat(" - Use simple random sampling to select 20 towns per country, each town equally likely.\n")
cat(" - Collect drinking-water samples in the same way at each selected town, then assay calcium (ppm).\n\n")

cat("Q2. Compare calcium means: Canada vs USA\n")
cat("Hypotheses: H0: mu_CAN - mu_USA = 0 ; H1: not equal\n")
cat(sprintf("Mean(Canada) = %.3f ppm; SD = %.3f; n = %d\n", mean(x_can), sd(x_can), length(x_can)))
cat(sprintf("Mean(USA)    = %.3f ppm; SD = %.3f; n = %d\n", mean(y_usa), sd(y_usa), length(y_usa)))
cat(sprintf("Mean difference (CAN - USA) = %.3f ppm\n", mean_diff))
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f; 95%% CI: [%.3f, %.3f]\n",
            unname(t_welch$statistic), unname(t_welch$parameter), t_welch$p.value,
            t_welch$conf.int[1], t_welch$conf.int[2]))
if (!is.null(shap_can) && !is.null(shap_usa)) {
  cat(sprintf("Normality checks (Shapiro-Wilk): Canada p = %.3f, USA p = %.3f (p>0.05 supports approximate normality).\n",
              shap_can$p.value, shap_usa$p.value))
}
cat("Conclusion: At alpha = 0.05, ")
if (t_welch$p.value < 0.05) {
  cat("reject H0, there is evidence of a difference in mean calcium between Canada and the USA.\n\n")
} else {
  cat("fail to reject H0, insufficient evidence of a difference in mean calcium between Canada and the USA.\n\n")
}

cat("Q3. Mortality vs Calcium (n = 25 towns)\n")
cat("(a) Scatter plot saved to results/figures/scatter_mortality_calcium.png\n")
cat(sprintf("(b) Pearson correlation: r = %.3f, p = %.4f (H0: rho = 0)\n",
            unname(cor_pearson$estimate), cor_pearson$p.value))
cat("Interpretation: p<0.05 indicates evidence of linear association; sign of r shows direction.\n")
cat("(c) Linear regression of mortality on calcium:\n")
s <- summary(fit)
cat(sprintf("  mortality = %.3f + (%.3f) * calcium\n",
            coef(fit)[1], coef(fit)[2]))
cat(sprintf("  Slope p-value = %.4f; R-squared = %.3f\n",
            s$coefficients["calcium","Pr(>|t|)"], s$r.squared))
cat(sprintf("  Prediction at calcium = 40 ppm (95%% PI): fit = %.2f, lower = %.2f, upper = %.2f\n",
            pred40[1], pred40[2], pred40[3]))
cat(sprintf("  Reliability note: 40 ppm is %s the observed calcium range [%.1f, %.1f].\n",
            if (inside_range) "inside" else "outside", calcium_range[1], calcium_range[2]))
cat("  Predictions are more reliable for towns with calcium within the observed range and when residuals are roughly homoscedastic and normal.\n\n")

cat("Files written:\n")
cat(" - results/figures/boxplot_calcium_by_country.png\n")
cat(" - results/figures/scatter_mortality_calcium.png\n")
cat(" - results/summary.txt\n")
sink()
