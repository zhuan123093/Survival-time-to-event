# Generate Gaussian CRE/BRE, Gamma CRE, and stratified CAR data.
generate_data <- function(trial, g) {
    n <- as.integer(g$n); car <- g$design == 'CAR'
    seed <- if (car) as.integer(90517L + 3000000L + ifelse(g$allocation == .5, 0L, 20000000L) + trial) else
        as.integer(41027000L + c(CRE = 0L, BRE = 80000000L)[g$design] +
            ifelse(g$allocation == .5, 0L, 40000000L) + n * 20000L + trial +
            ifelse(g$experiment == 'gamma', 900000000L, 0L))
    set.seed(seed)
    if (car) {
        biological_X <- matrix(rnorm(n * 30L), n)
        ss <- factor(rep(1:3, each = 50L), levels = 1:3)
        biological_X[, 2] <- qnorm((as.integer(ss) - 1 + runif(n))/3)
        A <- integer(n)
        for (s in 1:3) {
            ii <- which(ss == s); n1 <- ceiling(50 * g$allocation)
            A[ii] <- sample(c(rep.int(1L, n1), rep.int(0L, 50L - n1)))
        }
        entry <- runif(n, 0, 6); base_exp <- rexp(n)
        signal <- c(-.4, .5, .2)[as.integer(ss)] * rowSums(biological_X[, 1:10, drop = FALSE])
        set.seed(seed + 10000000L); W <- matrix(rnorm(n * 30L), n)
        X <- cbind(biological_X[, 1:10, drop = FALSE], W[, 1:5, drop = FALSE])
    } else {
        X <- if (g$experiment == 'gamma') matrix((rgamma(n * 25L, shape = 2, scale = 1) - 2)/sqrt(2), n) else
            matrix(rnorm(n * 25L), n)
        A <- if (g$design == 'CRE') sample(c(rep.int(1L, ceiling(n * g$allocation)),
            rep.int(0L, n - ceiling(n * g$allocation)))) else rbinom(n, 1L, g$allocation)
        ss <- factor(rep(1L, n))
        entry <- runif(n, 0, 6); base_exp <- rexp(n)
        signal <- .4 * rowSums(X[, 1:15, drop = FALSE])
    }
    colnames(X) <- paste0('x', seq_len(ncol(X)))
    event_time <- base_exp/((log(2)/12) * exp(signal + g$logHR * A))
    observed <- administrative_censoring(entry, event_time, list(event_fraction = .7))
    c(list(X = X, A = A, stratum = ss, entry = entry, base_exp = base_exp,
        signal = signal, seed = seed), observed)
}

compute_trial <- function(trial, g) {
    d <- generate_data(trial, g); ps <- as.integer(strsplit(g$ps, ',', fixed = TRUE)[[1]])
    n <- length(d$A); car <- g$design == 'CAR'
    der <- if (car) compute_stratified_derived_outcomes(d$time, d$status, d$A, d$stratum) else
        compute_derived_outcomes(d$time, d$status, d$A)
    lr <- if (car) c(u = der$score/n, v = der$information/n) else fit_logrank(d$time, d$status, d$A)[c('u', 'v')]
    pooled <- pooled_grid(d, ps); loo <- chiang_fit(d, ps)
    rows <- list()
    for (j in seq_along(ps)) {
        robin <- shared_robin(d$X[, seq_len(ps[j]), drop = FALSE], d$A, der$O0, der$O1, der,
            if (car) d$stratum else NULL)
        fits <- list(lr = lr, ols = robin[c('u_ols', 'v_ols')], dof = robin[c('u_dof', 'v_dof')],
            loo = loo$parts[j, c('u_check', 'v')], cloo = pooled$parts[j, 'cloo', ],
            ploo = pooled$parts[j, 'ploo', ], oracle = pooled$parts[j, 'oracle', ])
        rows[[j]] <- unlist(lapply(names(fits), function(method) {
            uv <- unname(fits[[method]]); z <- if (all(is.finite(uv)) && uv[2] > 0) sqrt(n) * uv[1]/sqrt(uv[2]) else NA_real_
            setNames(c(uv, z), paste0(c('u_', 'v_', 'z_'), method))
        }))
    }
    raw <- do.call(rbind, rows); rownames(raw) <- as.character(ps)
    stopifnot(sum(d$status) == ceiling(.7 * n), d$entry_gap > 0)
    list(raw = raw, parts = loo$parts, between = pooled$between,
        meta = c(seed = d$seed, n1 = sum(d$A), events = sum(d$status), cutoff = d$cutoff))
}
