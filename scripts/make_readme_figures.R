# Draws figures/fig0_records_vs_breeding_values.png, fig1_ghs_genetic_variance.png and fig2_where_variance_goes.png
# for README.md from the pooled
# step-1 results in results/ (SIMULATED data). Run from the repository root:  Rscript scripts/make_readme_figures.R
# Needs ggplot2 and ragg. Colours: categorical slots checked for colour-vision deficiency (adjacent pairs).
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[1] else "."
suppressPackageStartupMessages({ library(ggplot2); library(ragg) })

surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
c_def <- "#2a78d6"; c_pc <- "#eb6834"                       # categorical slots 1 and 2
c_g <- "#2a78d6"; c_h <- "#eb6834"; c_s <- "#1baf7a"; c_e <- "#eda100"   # slots 1-4, adjacent order validated

theme_step1 <- theme_minimal(base_size = 15) +
  theme(plot.background = element_rect(fill = surface, colour = NA),
        panel.background = element_rect(fill = surface, colour = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = grid, linewidth = 0.4),
        text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(colour = ink2, size = 13),
        plot.caption = element_text(colour = ink2, size = 11, hjust = 0),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.margin = margin(14, 18, 10, 14))

# ---- Figure 1: GHS sigma2_g per data set, default prior vs PC prior -------------------------------
d0 <- read.csv(file.path(repo, "results/step1_default_prior/sb_by_rep.csv"))
d1 <- read.csv(file.path(repo, "results/step1_pc_prior/sb_by_rep.csv"))
pick <- function(d, lab) { x <- d[d$model == "GHS" & !is.na(d$s2g), c("ncow", "s2g", "collapsed")]; x$prior <- lab; x }
d <- rbind(pick(d0, "INLA's default prior"), pick(d1, "PC prior"))
d$prior <- factor(d$prior, levels = c("INLA's default prior", "PC prior"))
d$grp <- factor(sprintf("%s cows", format(d$ncow, big.mark = ",")), levels = c("1,000 cows", "2,000 cows", "4,000 cows"))
d$xpos <- as.numeric(d$grp) + ifelse(d$prior == "PC prior", 0.2, -0.2)
d$xpos[d$ncow == 4000] <- 3                                   # no PC refit at 4,000 cows
set.seed(7); d$xj <- d$xpos + runif(nrow(d), -0.13, 0.13)
m <- aggregate(s2g ~ prior + grp + xpos, data = d, FUN = mean)
m$n <- aggregate(s2g ~ prior + grp + xpos, data = d, FUN = length)$s2g
m$bias <- 100 * (m$s2g / 0.06 - 1)
pct <- function(b) ifelse(abs(b) < 0.5, "0 %", sub("-", "−", sprintf("%+.0f %%", b)))
m$lab <- sprintf("mean %.3f\n%s", m$s2g, pct(m$bias))
cat("Figure 1 means:\n"); print(m)
ncol0 <- sum(d$collapsed & d$prior == "INLA's default prior" & d$ncow == 1000)

p1 <- ggplot(d, aes(xj, s2g)) +
  geom_hline(yintercept = 0.06, linetype = "22", colour = ink, linewidth = 0.6) +
  annotate("text", x = 2.64, y = 0.0625, label = "truth 0.060", vjust = 0, colour = ink, size = 4.3) +
  geom_point(aes(fill = prior), shape = 21, size = 2.9, colour = surface, stroke = 0.5, alpha = 0.9) +
  geom_segment(data = m, aes(x = xpos - 0.17, xend = xpos + 0.17, y = s2g, yend = s2g), colour = ink, linewidth = 1.1,
               lineend = "round") +
  geom_text(data = m, aes(x = xpos, y = 0.142, label = lab), colour = ink, size = 3.8, lineheight = 0.95, vjust = 0) +
  annotate("text", x = 0.8, y = -0.009, label = sprintf("%d of 40 collapsed to about 0", ncol0), colour = ink2,
           size = 3.8) +
  scale_fill_manual(values = c(c_def, c_pc)) +
  scale_x_continuous(breaks = 1:3, labels = levels(d$grp), limits = c(0.45, 3.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.14, 0.02), limits = c(-0.014, 0.162)) +
  labs(x = NULL, y = expression(hat(sigma)[g]^2),
       title = "Genetic variance estimated by GHS, the true model: one dot per data set",
       subtitle = "Step 1, one record per cow. 10 replicates x 4 herd structures at 1,000 and 2,000 cows; 5 data sets at 4,000.",
       caption = "SIMULATED data. The 4,000-cow data sets were not refitted with the PC prior. Bar = mean; percentage = bias against 0.060.") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_step1
agg_png(file.path(repo, "figures/fig1_ghs_genetic_variance.png"), width = 1800, height = 1050, res = 150, background = surface)
print(p1); invisible(dev.off())

# ---- Figure 2: where each model puts the variance, 2,000 cows ------------------------------------
v <- d0[d0$ncow == 2000, ]
mm <- function(mod, col) { x <- v[v$model == mod, col]; if (all(is.na(x))) 0 else mean(x, na.rm = TRUE) }
models <- c("G", "GH", "GS", "GHS")
comp <- rbind(data.frame(model = "Truth", term = c("genetic", "herd", "field", "residual"), value = c(0.06, 0.15, 0.41, 0.29)),
              do.call(rbind, lapply(models, function(mo) data.frame(model = mo, term = c("genetic", "herd", "field", "residual"),
                value = c(mm(mo, "s2g"), mm(mo, "s2h"), mm(mo, "s2s"), mm(mo, "s2e"))))))
comp$value[is.na(comp$value)] <- 0
cat("\nFigure 2 components (2,000 cows):\n"); print(reshape(comp, idvar = "model", timevar = "term", direction = "wide"))
comp$model <- factor(comp$model, levels = rev(c("Truth", models)))
comp$term <- factor(comp$term, levels = c("residual", "field", "herd", "genetic"))   # genetic drawn first (left)
comp <- comp[comp$value > 0, ]
comp <- comp[order(comp$model, -as.numeric(comp$term)), ]
comp$end <- ave(comp$value, comp$model, FUN = cumsum); comp$mid <- comp$end - comp$value / 2
comp$lab <- ifelse(comp$value >= 0.1, sprintf("%s %.2f", comp$term, comp$value), "")
gl <- comp[comp$term == "genetic", ]
gl$tot <- tapply(comp$value, comp$model, sum)[as.character(gl$model)]
gl$txt <- ifelse(gl$model == "Truth", sprintf("genetic %.3f", gl$value),
                 sprintf("genetic %.3f (%s)", gl$value, sub("-", "−", sprintf("%+.0f %%", 100 * (gl$value / 0.06 - 1)))))
xmax <- max(gl$tot) + 0.36

p2 <- ggplot(comp, aes(y = model)) +
  geom_col(aes(x = value, fill = term), colour = surface, linewidth = 0.9, width = 0.62, orientation = "y") +
  geom_text(aes(x = mid, label = lab), colour = ink, size = 3.8) +
  geom_text(data = gl, aes(x = tot + 0.02, label = txt), hjust = 0, colour = ink, size = 3.9) +
  scale_fill_manual(values = c(genetic = c_g, herd = c_h, field = c_s, residual = c_e),
                    breaks = c("genetic", "herd", "field", "residual")) +
  scale_x_continuous(limits = c(0, xmax), breaks = seq(0, 1.2, 0.2), expand = c(0, 0)) +
  labs(x = "variance (mean of 40 data sets)", y = NULL,
       title = "Where each model puts the variance of a term it leaves out (2,000 cows)",
       subtitle = "G and GS put the missing herd and field into the residual; GH puts the field into the herd. None puts it into genetics.",
       caption = paste("SIMULATED data, INLA's default prior. Field = SPDE variance parameter (truth 0.41; about 0.35 across the herds).",
                       "Genetic variance is low in every model because of the prior (see the PC refit).", sep = "\n")) +
  theme_step1 + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid, linewidth = 0.4),
                      axis.text.y = element_text(colour = ink, size = 14, face = "bold"))
agg_png(file.path(repo, "figures/fig2_where_variance_goes.png"), width = 1800, height = 900, res = 150, background = surface)
print(p2); invisible(dev.off())
# ---- Figure 0: predicting records is not predicting breeding values (cross-validation, 2,000 cows) --------
# Source: results/step1_default_prior/sb_cv_testpop_by_rep.csv (one row per data set x model; each measure computed
# within each of the 5 test groups and averaged, as in README section 4). Only data sets with all 5 groups are used.
cv <- read.csv(file.path(repo, "results/step1_default_prior/sb_cv_testpop_by_rep.csv"))
cv <- cv[cv$n_folds == 5 & grepl("^n2000", cv$cell), ]
f0 <- aggregate(cbind(r_ebv_a, cor_y_yhat) ~ model, cv, mean)
f0$n <- aggregate(r_ebv_a ~ model, cv, length)$r_ebv_a
cat("\nFigure 0 (2,000 cows, test cows):\n"); print(f0)
lab_r <- "Accuracy of breeding values, r(EBV, a)"
lab_y <- "Prediction of records, cor(y, ŷ)"
long <- rbind(data.frame(model = f0$model, what = lab_r, v = f0$r_ebv_a),
              data.frame(model = f0$model, what = lab_y, v = f0$cor_y_yhat))
long$what <- factor(long$what, levels = c(lab_r, lab_y))
g_r <- f0$r_ebv_a[f0$model == "GHS"] - f0$r_ebv_a[f0$model == "G"]
g_y <- f0$cor_y_yhat[f0$model == "GHS"] - f0$cor_y_yhat[f0$model == "G"]
lev <- rev(c("G", "GH", "GS", "GHS"))
f0$model <- factor(f0$model, levels = lev); long$model <- factor(long$model, levels = lev)
ends <- long[long$model %in% c("G", "GHS"), ]
ends$hj <- ifelse(ends$what == lab_r, 1.45, -0.45)
p0 <- ggplot(long, aes(y = model)) +
  geom_segment(data = f0, aes(x = r_ebv_a, xend = cor_y_yhat, yend = model), colour = grid, linewidth = 2.2,
               lineend = "round") +
  geom_point(aes(x = v, fill = what), shape = 21, size = 5, colour = surface, stroke = 1) +
  geom_text(data = ends, aes(x = v, label = sprintf("%.2f", v), hjust = hj), colour = ink, size = 4.2) +
  scale_fill_manual(values = c(c_def, c_pc)) +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
  labs(x = "correlation, for cows predicted without their own records", y = NULL,
       title = "Herd and field raise record prediction five times more than breeding-value accuracy",
       subtitle = sprintf("From G to GHS: prediction of records +%.2f, accuracy of breeding values +%.2f. Step 1, 2,000 cows, five-fold cross-validation.",
                          g_y, g_r),
       caption = sprintf(paste0("SIMULATED data, INLA's default prior. Each dot: mean of %d data sets (replicates 1-5 x 4 herd structures); ",
                                "measures computed within each test group.\nG = genetic only; GH + herd; GS + spatial field; ",
                                "GHS herd + field (the true model)."), min(f0$n))) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_step1 + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(colour = grid, linewidth = 0.4),
                      axis.text.y = element_text(colour = ink, size = 14, face = "bold"))
agg_png(file.path(repo, "figures/fig0_records_vs_breeding_values.png"), width = 1800, height = 760, res = 150,
        background = surface)
print(p0); invisible(dev.off())
cat("written\n")
