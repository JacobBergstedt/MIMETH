
glmnet_logreg <- function(x, y, q, ...) {
  fit <- suppressWarnings(glmnet(x, y, family = "binomial", pmax = q, ...))
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  list(selected = ret, path = sequence)
}


cv_stability_selection <- function (cv_folds, meth, y, alpha, q, tol) {
  
  
  map_dfr(1:10, stability_selection_fold, cv_folds, meth, y, alpha = alpha, q = q, tol = tol) %>%
    bind_rows()
  
  
}

stability_selection_fold <- function(fold, cv_folds, meth, y, alpha, q, tol) {
  
  id_train <- cv_folds != fold
  id_test <- cv_folds == fold
  
  meth_train <- meth[id_train, ]
  y_train <- y[id_train]
  
  meth_test <- meth[id_test, ]
  y_test <- y[id_test]
  
  selection <- stabsel(meth_train, 
                       y_train,
                       fitfun = glmnet_logreg,
                       args.fitfun = list(intercept = TRUE, alpha = alpha),
                       q = q,
                       PFER = tol,
                       papply = mclapply,
                       mc.cores = 10)
  
  df <- as.data.frame(meth_train[, selection$selected])
  df$y <- y_train
  m <- glm(y ~ ., df, family = "binomial")
  
  meth_test <- as.data.frame(meth_test[, selection$selected])
  y_pred <- predict(m, meth_test, type = "response")
  
  class_pred <- ifelse(y_pred > 0.50, "Yes", "No")
  class_pred <- factor(class_pred, c("No", "Yes"))
  levels(class_pred) <- levels(y_test)
  
  tibble(pheno = "CMV_serostatus", 
         classrate = mean(class_pred == y_test),
         nr_predictors = length(selection$selected),
         Fold = fold)
  
}


est_pred_accuracy <- function(y_test, y_train, meth_test, meth_train) {
  
  df <- as.data.frame(meth_train)
  df$y <- y_train
  m <- glm(y ~ ., df, family = "binomial")
  y_pred <- predict(m, as.data.frame(meth_test), type = "response")
  
  class_pred <- ifelse(y_pred > 0.50, "Yes", "No")
  class_pred <- factor(class_pred, c("No", "Yes"))
  levels(class_pred) <- levels(y_test)
  
  tibble(pheno = "CMV_serostatus", 
         classrate = sum(class_pred == y_test) / length(y_test), 
         nr_predictors = ncol(meth_test))
  
  
}


cv_alpha_logreg <- function(X, y, alpha_seq, fold_list, n_lambda = 1e2, n_cores) {
  
  map(alpha_seq, function(alpha) cv_logreg(X, y, alpha = alpha, fold_list = fold_list, n_lambda = n_lambda))
}

fit_glmnet_for_fold <- function(fold, cv_folds, X, y, alpha, lambda) {
  
  
  ids <- fold != cv_folds
  
  X_train <- X[ids, ]
  y_train <- y[ids]
  
  X_test <- X[!ids, ]
  y_test <- y[!ids]
  
  fit <- glmnet(X_train, y_train, family = "binomial", alpha = alpha, lambda = lambda)
  pred_class <- predict(fit, newx = X_test, type = "class")
  pred_prob <- predict(fit, newx = X_test, type = "response")
  pred_nonzero <- predict(fit, type = "nonzero")
  CR <- colMeans(y_test == pred_class)
  
  list(alpha = alpha,
       lambda = lambda,
       lambda_nr = seq_along(lambda),
       fold_nr = cv_folds,
       observed_class = y_test,
       CR = CR,
       pred_class = pred_class,
       pred_prob = pred_prob,
       pred_nonzero = map_int(pred_nonzero, length))
  
}



cv_logreg <- function(X, y, alpha, fold_list, n_lambda) {
  
  cv_folds_logreg <- function(cv_folds, X, y, lambda) {
    
    mclapply(unique(cv_folds), function(fold) fit_glmnet_for_fold(fold, cv_folds = cv_folds, X = X, y = y, alpha = alpha, lambda = lambda), mc.cores = 10)
    
  }
  
  lambda <- glmnet(X, y, alpha = alpha, family = "binomial", nlambda = n_lambda)$lambda
  map(fold_list, cv_folds_logreg, X = X, y = y, lambda = lambda)
  
}

plot_cv <- function(res) {
  res %>% 
    ggplot(aes(x = Lambda_nr, y = Acc)) +
    geom_boxplot() +
    facet_wrap(vars(Alpha))
  
}