# OLS/dof: fitted slopes and finite-sample multiplier.

krotka_gamma <- function(n1, n0, k) {
    n <- n1 + n0
    if (n1 <= k + 2 || n0 <= k + 2) 
        return(NA_real_)
    df_part <- n/(n - (n0/n1 + n1/n0 - 1) * k - (n0/n1 + n1/n0))
    vif_part <- 1 + k * ((n0/n)^2/(n1 - k - 2) + (n1/n)^2/(n0 - k - 2))
    df_part * vif_part
}

# Group-specific OLS adjustment and degrees-of-freedom correction.

shared_robin <- function(X, A, O0, O1, der, stratum = NULL) {
    n <- length(A)
    p <- ncol(X)
    empty <- c(u_ols = NA_real_, v_ols = NA_real_, u_dof = NA_real_, v_dof = NA_real_, gamma = NA_real_, 
        fit_failed = 1)
    tryCatch({
        beta <- vector("list", 2)
        for (a in 0:1) {
            idx <- which(A == a)
            x <- X[idx, , drop = FALSE]
            y <- if (a == 1) 
                O1[idx]
            else O0[idx]
            if (is.null(stratum)) {
                x <- scale(x, center = TRUE, scale = FALSE)
                beta[[a + 1L]] <- lm.fit(x, y, singular.ok = FALSE)$coefficients
            }
            else {
                cross_x <- cross_y <- list()
                ss <- as.integer(stratum[idx])
                for (s in unique(ss)) {
                  ii <- which(ss == s)
                  xx <- scale(x[ii, , drop = FALSE], center = TRUE, scale = FALSE)
                  yy <- y[ii] - mean(y[ii])
                  cross_x[[length(cross_x) + 1L]] <- crossprod(xx)
                  cross_y[[length(cross_y) + 1L]] <- crossprod(xx, yy)
                }
                beta[[a + 1L]] <- solve(Reduce("+", cross_x), Reduce("+", cross_y))
            }
        }
        if (is.null(stratum)) {
            xall <- rbind(X[A == 0L, , drop = FALSE], X[A == 1L, , drop = FALSE])
            bar <- colMeans(xall)
            xc <- sweep(X, 2L, bar, "-")
            covx <- cov(xall)
        }
        else {
            xc <- X
            covs <- list()
            for (s in levels(stratum)) {
                ii <- which(stratum == s)
                xall <- rbind(X[ii[A[ii] == 0L], , drop = FALSE], X[ii[A[ii] == 1L], , drop = FALSE])
                xc[ii, ] <- sweep(X[ii, , drop = FALSE], 2L, colMeans(xall), "-")
                covs[[s]] <- length(ii)/n * cov(xall)
            }
            covx <- Reduce("+", covs)
        }
        u <- der$score/n - (sum(xc[A == 1L, , drop = FALSE] %*% beta[[2L]]) - sum(xc[A == 0L, , drop = FALSE] %*% 
            beta[[1L]]))/n
        bs <- beta[[1L]] + beta[[2L]]
        q <- mean(A)
        v <- der$information/n - q * (1 - q) * as.numeric(t(bs) %*% covx %*% bs)
        gamma <- krotka_gamma(sum(A), sum(1 - A), p)
        c(u_ols = u, v_ols = v, u_dof = u, v_dof = gamma * v, gamma = gamma, fit_failed = 0)
    }, error = function(e) empty)
}

# Moore-Penrose convention for the centered pooled design.

mf_inverse <- function(X) {
    G <- crossprod(X)
    cc <- rcond(G)
    p <- ncol(X)
    if (is.finite(cc) && cc > 1e-09) 
        return(list(W = chol2inv(chol(G)), rank = p))
    f <- svd(X)
    keep <- f$d > max(dim(X)) * .Machine$double.eps * max(f$d)
    W <- if (any(keep)) 
        tcrossprod(sweep(f$v[, keep, drop = FALSE], 2, f$d[keep], `/`))
    else matrix(0, p, p)
    list(W = W, rank = sum(keep))
}

# Known conditional variances for the three Gaussian tertiles.

oracle_CAR_moments <- function() {
    boundaries <- qnorm(c(0, 1/3, 2/3, 1))
    rows <- list()
    for (k in 1:3) {
        lower <- boundaries[k]
        upper <- boundaries[k + 1]
        prob <- 1/3
        mean <- (dnorm(lower) - dnorm(upper))/prob
        left <- if (is.finite(lower)) 
            lower * dnorm(lower)
        else 0
        right <- if (is.finite(upper)) 
            upper * dnorm(upper)
        else 0
        variance <- 1 + (left - right)/prob - mean^2
        rows[[k]] <- data.frame(stratum = k, lower = lower, upper = upper, mean = mean, variance = variance)
    }
    do.call(rbind, rows)
}

# Symmetrized pair-kernel square evaluated by a low-rank identity.

bd_pair <- function(F, G, dense = FALSE) {
    m <- nrow(F)
    if (dense) {
        K <- tcrossprod(F, G)
        diag(K) <- 0
        return(sum((K + t(K))^2)/(2 * m))
    }
    gf <- crossprod(G, F)
    dg <- rowSums(F * G)
    z <- (sum(crossprod(F) * crossprod(G)) + sum(gf * t(gf)) - 2 * sum(dg^2))/m
    if (!is.finite(z) || z < 0) {
        K <- tcrossprod(F, G)
        diag(K) <- 0
        z <- sum((K + t(K))^2)/(2 * m)
    }
    z
}

# Ordered off-diagonal sum of (O_i1+O_i0) C_ij C_ji (O_j1+O_j0), C = Z X'.
off_diagonal_product <- function(Z, X, response) {
    diagonal <- rowSums(Z * X)
    moment <- crossprod(X, Z * response)
    sum(moment * t(moment)) - sum((diagonal * response)^2)
}

# Components I and II follow the displayed PLOO/CLOO/oracle formulas.
# Gram matrices and response moments both use divisor n; their scales cancel
# in the fitted sample weights. Known-Sigma weights retain the divisor n.
pooled_fit <- function(X, A, O0, O1, vlr, stratum = 1L, CAR = FALSE,
                       methods = c('cloo', 'ploo', 'oracle')) {
    m <- nrow(X); p <- ncol(X); q <- mean(A); v <- q * (1 - q)
    pa <- ifelse(A == 1, q, 1 - q)
    yy <- ifelse(A == 1, O1, O0); sig <- O0 + O1
    rows <- list()
    if (any(c('cloo', 'ploo') %in% methods)) {
        inverse <- mf_inverse(X); U <- X %*% inverse$W
        h <- rowSums(U * X)
        C <- U
        safe <- inverse$rank == p & (1 - h) > 1e-6
        C[safe, ] <- sweep(C[safe, , drop = FALSE], 1, 1 - h[safe], '/')
        for (i in which(!safe)) C[i, ] <- drop(X[i, , drop = FALSE] %*% mf_inverse(X[-i, , drop = FALSE])$W)
        rows <- list(cloo = C, ploo = U)
    }
    sigma <- rep(1, p)
    if (CAR && p >= 2L) sigma[2] <- oracle_CAR_moments()$variance[stratum]
    oracle <- sweep(X, 2, m * sigma, '/')
    rows$oracle <- oracle
    rows <- rows[methods]
    score <- drop(crossprod(X, yy/pa))
    out <- matrix(0, length(rows), 4, dimnames = list(names(rows), c('u', 'I', 'II', 'v')))
    for (method in names(rows)) {
        Z <- rows[[method]]; hz <- rowSums(Z * X)
        pred <- drop(Z %*% score) - hz * yy/pa
        phi <- A * O1 - (1 - A) * O0 + (q - A) * pred
        off_prediction <- drop(Z %*% crossprod(X, sig)) - hz * sig
        product_sum <- off_diagonal_product(Z, X, sig)
        if (method == 'ploo') {
            # I: (1+h_i+h_j+H_ij) H_ij times the two endpoint response sums.
            ordered_sum <- sum((1 + 2 * h) * sig * off_prediction) + product_sum
        } else if (method == 'cloo') {
            # C_ij = H_ij/(1-h_i) on full-rank deleted designs. The three
            # contributions are 2 C_ij, -sum_l C_li C_lj, and C_ij C_ji.
            column_squares <- rowSums((X %*% crossprod(Z)) * X) - hz^2
            triple_sum <- sum(off_prediction^2) - sum(sig^2 * column_squares)
            ordered_sum <- 2 * sum(sig * off_prediction) - triple_sum + product_sum
        } else {
            # I: (1+H_ij^oracle) H_ij^oracle, omitting the 2/n correction.
            ordered_sum <- sum(sig * off_prediction) + product_sum
        }
        I <- vlr - v/m * ordered_sum
        # II: arm-weighted pair squares. For CLOO, the symmetric pair is
        # H_ij {O_i/(1-h_j) + O_j/(1-h_i)}. bd_pair sums each pair once / n.
        II <- (1 - q) * bd_pair(Z, X * O1) + q * bd_pair(Z, X * O0)
        out[method, ] <- c(mean(phi), I, II, I + II)
    }
    out
}

pooled_grid <- function(d, ps, methods = c('cloo', 'ploo', 'oracle')) {
    n <- length(d$A); ans <- array(0, c(length(ps), length(methods), 2),
        dimnames = list(as.character(ps), methods, c('u', 'v')))
    components <- array(0, c(length(ps), length(methods), 3),
        dimnames = list(as.character(ps), methods, c('I', 'II', 'between')))
    means <- weights <- numeric()
    for (k in unique(as.integer(d$stratum))) {
        ii <- which(as.integer(d$stratum) == k); m <- length(ii); A <- d$A[ii]; q <- mean(A)
        der <- compute_derived_outcomes(d$time[ii], d$status[ii], A)
        means <- c(means, q * mean(der$O1[A == 1]) - (1 - q) * mean(der$O0[A == 0]))
        weights <- c(weights, m/n)
        Xc <- scale(d$X[ii, seq_len(max(ps)), drop = FALSE], TRUE, FALSE)
        for (j in seq_along(ps)) {
            fit <- pooled_fit(Xc[, seq_len(ps[j]), drop = FALSE], A, der$O0, der$O1,
                der$information/m, k, length(unique(d$stratum)) > 1L, methods)
            ans[j, , ] <- ans[j, , ] + m/n * fit[, c('u', 'v')]
            components[j, , c('I', 'II')] <- components[j, , c('I', 'II')] + m/n * fit[, c('I', 'II')]
        }
    }
    between <- sum(weights * (means - sum(weights * means))^2)
    ans[, , 'v'] <- ans[, , 'v'] + between
    components[, , 'between'] <- between
    list(parts = ans, components = components, between = between)
}

