# Stage 04: fit GAM + lasso point-process models per training set; spatial block CV;
# permutation importance; response curves
source("pipeline/00_config.R"); source("pipeline/functions.R")
md <- readRDS(file.path(CFG$data_dir, "model_data.rds"))
dir.create("models", showWarnings = FALSE)
vars <- md$vars

cv_rows <- list(); imp_rows <- list(); rc_rows <- list()
for (s in names(md$sets)) {
  message("== ", s)
  S <- md$sets[[s]]; d <- S$d

  # spatial block cross-validation
  for (k in sort(unique(d$fold))) {
    tr <- d[d$fold != k, ]; te <- d[d$fold == k, ]
    if (sum(te$pa == 1) < 5) next
    m <- fit_set(tr, vars, S$bias, CFG)
    pi <- predict_index(m, te)
    for (a in c("gam", "lasso", "ens")) {
      p <- pi[[a]][te$pa == 1]; b <- pi[[a]][te$pa == 0]
      cv_rows[[length(cv_rows) + 1]] <- data.table(
        set = s, algo = a, fold = k, n_test_pres = length(p),
        auc = auc(p, b), boyce = boyce(p, b),
        test_pres_above_p10 = mean(p >= 0.10))
    }
  }

  # final model on all data
  m <- fit_set(d, vars, S$bias, CFG)
  saveRDS(m, file.path("models", paste0(s, ".rds")))

  # permutation importance (drop in ensemble AUC)
  set.seed(CFG$seed)
  ev <- rbind(d[d$pa == 1, ], d[d$pa == 0, ][sample(sum(d$pa == 0), min(10000, sum(d$pa == 0))), ])
  base_auc <- auc(predict_index(m, ev)$ens[ev$pa == 1], predict_index(m, ev)$ens[ev$pa == 0])
  for (v in vars) {
    drops <- vapply(1:3, function(r) {
      e2 <- ev; e2[[v]] <- sample(e2[[v]])
      ip <- predict_index(m, e2)$ens
      base_auc - auc(ip[e2$pa == 1], ip[e2$pa == 0])
    }, 1)
    imp_rows[[length(imp_rows) + 1]] <- data.table(set = s, variable = v, auc_drop = mean(drops))
  }

  # response curves: other predictors at presence medians, effort fixed at reference
  med <- as.data.frame(lapply(d[d$pa == 1, c(vars, S$bias), drop = FALSE], median))
  for (v in vars) {
    xs <- quantile(d[[v]][d$pa == 0], c(0.01, 0.99), names = FALSE)
    nd <- med[rep(1, 100), , drop = FALSE]; nd[[v]] <- seq(xs[1], xs[2], length.out = 100)
    for (a in c("gam", "lasso")) {
      eta <- predict_eta(a, m[[a]], nd, m$fix_bias)
      rc_rows[[length(rc_rows) + 1]] <- data.table(set = s, variable = v, algo = a,
                                                   x = nd[[v]], relative = exp(eta - max(eta)))
    }
  }
}

cv <- rbindlist(cv_rows)
fwrite(cv, file.path(CFG$out_dir, "cv_folds.csv"))
cv_sum <- cv[, .(folds = .N, auc = mean(auc), auc_sd = sd(auc), boyce = mean(boyce, na.rm = TRUE),
                 boyce_sd = sd(boyce, na.rm = TRUE), test_pres_above_p10 = mean(test_pres_above_p10)),
             by = .(set, algo)]
fwrite(cv_sum, file.path(CFG$out_dir, "cv_summary.csv"))
imp <- rbindlist(imp_rows)[, share := pmax(auc_drop, 0) / sum(pmax(auc_drop, 0)), by = set]
fwrite(imp, file.path(CFG$out_dir, "importance.csv"))
fwrite(rbindlist(rc_rows), file.path(CFG$out_dir, "response_curves.csv"))
print(cv_sum)
