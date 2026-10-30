# Current intercept CF + Chiang mHC3; HC3 is an internal term of mHC3.

cmo_variance <- function(A, r, P) {
    n <- length(A)
    n1 <- sum(A)
    n0 <- n - n1
    one <- A == 1
    zero <- !one
    h <- diag(P)
    K <- P^2
    diag(K) <- 0
    a <- r * A
    b <- r * (1 - A)
    base <- n/(n1 * (n1 - 1)) * sum(r[one]^2) + n/(n0 * (n0 - 1)) * sum(r[zero]^2)
    d1 <- (n1 - n0)/n1^2 * sum(h[one] * r[one]^2)
    d0 <- (n0 - n1)/n0^2 * sum(h[zero] * r[zero]^2)
    s11 <- n0^2 * n/n1^4 * sum(a * drop(K %*% a))
    s00 <- n1^2 * n/n0^4 * sum(b * drop(K %*% b))
    cross <- 2 * n/(n0 * n1) * sum(a * drop(K %*% b))
    c(HC3 = base, diagonal_1 = d1, diagonal_0 = d0, pair_11 = s11, pair_00 = s00, cross_10 = cross, mHC3 = base + 
        d1 + d0 + s11 + s00 - cross)
}

cmo_dense_grid <- function(X, A, O0, O1, ps, keep = FALSE) {
    n <- nrow(X)
    q <- mean(A)
    qa <- ifelse(A == 1, q, 1 - q)
    y <- ifelse(A == 1, O1, O0)
    arms <- lapply(0:1, function(a) {
        ii <- which(A == a)
        f <- qr(X[ii, , drop = FALSE], tol = max(dim(X)) * .Machine$double.eps)
        stopifnot(f$rank == ncol(X), identical(f$pivot, seq_len(ncol(X))))
        Q <- qr.Q(f, complete = TRUE)
        list(ii = ii, Q = Q, R = qr.R(f), score = drop(crossprod(Q, y[ii])))
    })
    fp <- qr(X, tol = max(dim(X)) * .Machine$double.eps)
    stopifnot(fp$rank == ncol(X), identical(fp$pivot, seq_len(ncol(X))))
    Qp <- qr.Q(fp)
    rows <- details <- list()
    for (p in ps) {
        r <- h <- numeric(n)
        pred <- matrix(0, n, 2)
        for (a in 0:1) {
            z <- arms[[a + 1]]
            rest <- seq.int(p + 1L, length(z$ii))
            N <- z$Q[, rest, drop = FALSE]
            complement <- rowSums(N^2)
            e <- drop(N %*% z$score[rest])
            r[z$ii] <- e/complement
            h[z$ii] <- 1 - complement
            beta <- backsolve(z$R[seq_len(p), seq_len(p), drop = FALSE], z$score[seq_len(p)])
            pred[, a + 1] <- drop(X[, seq_len(p), drop = FALSE] %*% beta)
            pred[z$ii, a + 1] <- y[z$ii] - r[z$ii]
        }
        P <- tcrossprod(Qp[, seq_len(p), drop = FALSE])
        vals <- cmo_variance(A, r * qa, P)
        u <- mean((2 * A - 1) * y - (A - q) * rowSums(pred))
        rows[[as.character(p)]] <- c(u_check = u, vals, max_arm_h = max(h), max_pool_h = max(diag(P)))
        if (keep) 
            details[[as.character(p)]] <- list(X = X[, seq_len(p), drop = FALSE], A = A, O0 = O0, O1 = O1, 
                residual = r, scaled_residual = r * qa, P = P, pred = pred, values = rows[[as.character(p)]])
    }
    list(parts = do.call(rbind, rows), details = details)
}

# Same mHC3 at large n, plus the saved saturated-arm NA convention.

cmo_fast_variance <- function(A, r, Q) {
    n <- length(A)
    n1 <- sum(A)
    n0 <- n - n1
    one <- A == 1
    zero <- !one
    h <- rowSums(Q^2)
    a <- r * A
    b <- r * (1 - A)
    M1 <- crossprod(Q, Q * a)
    M0 <- crossprod(Q, Q * b)
    base <- n/(n1 * (n1 - 1)) * sum(r[one]^2) + n/(n0 * (n0 - 1)) * sum(r[zero]^2)
    d1 <- (n1 - n0)/n1^2 * sum(h[one] * r[one]^2)
    d0 <- (n0 - n1)/n0^2 * sum(h[zero] * r[zero]^2)
    s11 <- n0^2 * n/n1^4 * (sum(M1^2) - sum(h^2 * a^2))
    s00 <- n1^2 * n/n0^4 * (sum(M0^2) - sum(h^2 * b^2))
    cross <- 2 * n/(n0 * n1) * sum(M1 * M0)
    c(HC3 = base, diagonal_1 = d1, diagonal_0 = d0, pair_11 = s11, pair_00 = s00, cross_10 = cross, mHC3 = base + 
        d1 + d0 + s11 + s00 - cross)
}

cmo_rank_case <- function(Z, A, O0, O1) {
    n <- nrow(Z)
    d <- ncol(Z)
    q <- mean(A)
    y <- ifelse(A == 1, O1, O0)
    pred <- matrix(0, n, 2)
    rank_bad <- FALSE
    maxh <- 0
    for (a in 0:1) {
        own <- which(A == a)
        X <- Z[own, , drop = FALSE]
        f <- qr(X, tol = max(dim(X)) * .Machine$double.eps)
        if (length(own) > d && f$rank == d) {
            Q <- qr.Q(f)
            h <- rowSums(Q^2)
            e <- qr.resid(f, y[own])
            beta <- qr.coef(f, y[own])
            pred[, a + 1] <- drop(Z %*% beta)
            pred[own, a + 1] <- y[own] - e/(1 - h)
            maxh <- max(maxh, h)
        }
        else {
            rank_bad <- TRUE
            maxh <- 1
            beta <- MASS::ginv(X, tol = max(dim(X)) * .Machine$double.eps) %*% y[own]
            pred[, a + 1] <- drop(Z %*% beta)
            for (i in own) {
                tr <- setdiff(own, i)
                beta <- MASS::ginv(Z[tr, , drop = FALSE], tol = max(dim(Z[tr, , drop = FALSE])) * .Machine$double.eps) %*% 
                  y[tr]
                pred[i, a + 1] <- sum(Z[i, ] * beta)
            }
        }
    }
    Q <- qr.Q(qr(Z))
    P <- tcrossprod(Q)
    r <- y - pred[cbind(seq_len(n), A + 1)]
    vals <- if (rank_bad) 
        setNames(rep(NA_real_, 7), c("HC3", "diagonal_1", "diagonal_0", "pair_11", "pair_00", "cross_10", 
            "mHC3"))
    else cmo_variance(A, r * ifelse(A == 1, q, 1 - q), P)
    values <- c(u_check = mean((2 * A - 1) * y - (A - q) * rowSums(pred)), vals, max_arm_h = maxh, max_pool_h = max(diag(P)))
    list(X = Z, A = A, O0 = O0, O1 = O1, residual = r, scaled_residual = r * ifelse(A == 1, q, 1 - q), 
        P = P, pred = pred, values = values, undefined_HC3 = rank_bad)
}

cmo_grid <- function(X, A, O0, O1, ps, keep = FALSE) {
    n <- nrow(X)
    ns <- c(sum(A == 0), sum(A == 1))
    full <- all(ns > ncol(X)) && all(vapply(0:1, function(a) qr(X[A == a, , drop = FALSE], tol = max(dim(X)) * 
        .Machine$double.eps)$rank == ncol(X), logical(1)))
    if (n <= 150 && full) 
        return(cmo_dense_grid(X, A, O0, O1, ps, keep))
    if (!full) {
        rows <- details <- list()
        for (p in ps) {
            Z <- X[, seq_len(p), drop = FALSE]
            safe <- all(ns > p) && all(vapply(0:1, function(a) qr(Z[A == a, , drop = FALSE], tol = max(dim(Z)) * 
                .Machine$double.eps)$rank == p, logical(1)))
            if (safe) {
                f <- cmo_dense_grid(Z, A, O0, O1, p, keep)
                rows[[as.character(p)]] <- f$parts[1, ]
                if (keep) 
                  details[[as.character(p)]] <- f$details[[1]]
            }
            else {
                f <- cmo_rank_case(Z, A, O0, O1)
                rows[[as.character(p)]] <- f$values
                if (keep) 
                  details[[as.character(p)]] <- f
            }
        }
        return(list(parts = do.call(rbind, rows), details = details))
    }
    q <- mean(A)
    qa <- ifelse(A == 1, q, 1 - q)
    y <- ifelse(A == 1, O1, O0)
    rows <- details <- list()
    for (p in ps) {
        Z <- X[, seq_len(p), drop = FALSE]
        r <- h <- numeric(n)
        pred <- matrix(0, n, 2)
        for (a in 0:1) {
            own <- which(A == a)
            f <- qr(Z[own, , drop = FALSE])
            Q <- qr.Q(f)
            hh <- rowSums(Q^2)
            e <- qr.resid(f, y[own])
            r[own] <- e/(1 - hh)
            h[own] <- hh
            beta <- qr.coef(f, y[own])
            pred[, a + 1] <- drop(Z %*% beta)
            pred[own, a + 1] <- y[own] - r[own]
        }
        Q <- qr.Q(qr(Z))
        vals <- cmo_fast_variance(A, r * qa, Q)
        values <- c(u_check = mean((2 * A - 1) * y - (A - q) * rowSums(pred)), vals, max_arm_h = max(h), 
            max_pool_h = max(rowSums(Q^2)))
        rows[[as.character(p)]] <- values
        if (keep) 
            details[[as.character(p)]] <- list(X = Z, A = A, O0 = O0, O1 = O1, residual = r, scaled_residual = r * 
                qa, P = tcrossprod(Q), pred = pred, values = values)
    }
    list(parts = do.call(rbind, rows), details = details)
}

# Combine current LOO within strata, retaining the original between term.
chiang_fit <- function (d, ps, keep = FALSE) 
{
    ans <- NULL
    details <- list()
    means <- weights <- numeric()
    for (k in unique(as.integer(d$stratum))) {
        ii <- which(as.integer(d$stratum) == k)
        m <- length(ii)
        n <- length(d$A)
        A <- d$A[ii]
        q <- mean(A)
        der <- compute_derived_outcomes(d$time[ii], d$status[ii], A)
        X <- cbind(1, scale(d$X[ii, seq_len(max(ps)), drop = FALSE], TRUE, FALSE))
        f <- cmo_grid(X, A, der$O0, der$O1, ps + 1L, keep)
        rownames(f$parts) <- as.character(ps)
        if (is.null(ans)) 
            ans <- 0 * f$parts
        cols <- setdiff(colnames(ans), c("max_arm_h", "max_pool_h"))
        ans[, cols] <- ans[, cols] + m/n * f$parts[, cols]
        ans[, "max_arm_h"] <- pmax(ans[, "max_arm_h"], f$parts[, "max_arm_h"])
        ans[, "max_pool_h"] <- pmax(ans[, "max_pool_h"], f$parts[, "max_pool_h"])
        means <- c(means, q * mean(der$O1[A == 1]) - (1 - q) * mean(der$O0[A == 0]))
        weights <- c(weights, m/n)
        if (keep) 
            for (j in seq_along(ps)) details[[paste(k, ps[j], sep = "_")]] <- f$details[[j]]
    }
    between <- sum(weights * (means - sum(weights * means))^2)
    list(parts = cbind(ans, between = between, v = ans[, "mHC3"] + between), details = details)
}
