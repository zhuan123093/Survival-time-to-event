# Current data utilities: event-driven censoring, derived outcomes and log-rank.

administrative_censoring <- function(entry, event_time, cfg) {
    n <- length(entry)
    calendar_event <- entry + event_time
    target <- if (is.null(cfg$event_fraction)) 
        NA_integer_
    else as.integer(ceiling(cfg$event_fraction * n))
    cutoff <- if (is.na(target)) 
        cfg$end_of_study
    else {
        if (target < 1L || target > n) 
            stop("Invalid target event count")
        sort(calendar_event, partial = target)[target]
    }
    gap <- cutoff - max(entry)
    if (!is.finite(gap) || gap <= 0) {
        stop(sprintf(paste0("Administrative cutoff %.6f precedes/completes at last ", "entry %.6f. Stop and ask before changing the design."), 
            cutoff, max(entry)))
    }
    status <- as.integer(calendar_event <= cutoff)
    if (!is.na(target) && sum(status) != target) 
        stop("Tied calendar events at cutoff")
    list(time = pmin(event_time, cutoff - entry), status = status, cutoff = cutoff, entry_gap = gap)
}

compute_derived_outcomes <- function(time, status, A) {
    n <- length(time)
    event_times <- sort(unique(time[status == 1L]))
    if (length(event_times) == 0L) {
        return(list(O0 = numeric(n), O1 = numeric(n), score = 0, information = 0, bridge_error = 0))
    }
    event_index <- match(time[status == 1L], event_times)
    d_all <- tabulate(event_index, nbins = length(event_times))
    d_1 <- tabulate(match(time[status == 1L & A == 1L], event_times), nbins = length(event_times))
    risk_total <- n - findInterval(event_times, sort(time), left.open = TRUE)
    time_1 <- sort(time[A == 1L])
    time_0 <- sort(time[A == 0L])
    risk_1 <- length(time_1) - findInterval(event_times, time_1, left.open = TRUE)
    risk_0 <- length(time_0) - findInterval(event_times, time_0, left.open = TRUE)
    if (any(risk_total <= 0L)) 
        stop("Empty risk set at an event time")
    w1 <- risk_0/risk_total
    w0 <- risk_1/risk_total
    breslow_increment <- d_all/risk_total
    cumulative_1 <- cumsum(w1 * breslow_increment)
    cumulative_0 <- cumsum(w0 * breslow_increment)
    last_event <- findInterval(time, event_times)
    compensator_1 <- ifelse(last_event > 0L, cumulative_1[pmax(last_event, 1L)], 0)
    compensator_0 <- ifelse(last_event > 0L, cumulative_0[pmax(last_event, 1L)], 0)
    failure_weight_1 <- numeric(n)
    failure_weight_0 <- numeric(n)
    is_event <- which(status == 1L)
    failure_weight_1[is_event] <- w1[match(time[is_event], event_times)]
    failure_weight_0[is_event] <- w0[match(time[is_event], event_times)]
    O1 <- failure_weight_1 - compensator_1
    O0 <- failure_weight_0 - compensator_0
    score <- sum(d_1 - d_all * risk_1/risk_total)
    information <- sum(d_all * risk_1 * risk_0/risk_total^2)
    score_from_O <- sum(A * O1 - (1L - A) * O0)
    list(O0 = O0, O1 = O1, score = score, information = information, bridge_error = score_from_O - score)
}

fit_logrank <- function(time, status, A, stratum = NULL) {
    dat <- data.frame(time = time, status = status, A = factor(A, levels = 0:1))
    strata <- survival::strata
    formula <- survival::Surv(time, status) ~ A
    if (!is.null(stratum)) {
        dat$stratum <- stratum
        formula <- survival::Surv(time, status) ~ A + strata(stratum)
    }
    fit <- survival::survdiff(formula, data = dat, rho = 0)
    score <- if (is.matrix(fit$obs)) 
        sum(fit$obs[2L, ] - fit$exp[2L, ])
    else fit$obs[2L] - fit$exp[2L]
    information <- fit$var[2L, 2L]
    n <- length(A)
    c(u = unname(score/n), v = unname(information/n), z = unname(score/sqrt(information)))
}

# Compute CAR risk sets within each stratum.

compute_stratified_derived_outcomes <- function(time, status, A, stratum) {
    n <- length(time)
    O0 <- numeric(n)
    O1 <- numeric(n)
    score <- 0
    information <- 0
    bridge_error <- 0
    for (level in levels(stratum)) {
        indices <- which(stratum == level)
        if (length(indices) == 0L) 
            next
        fit <- compute_derived_outcomes(time[indices], status[indices], A[indices])
        O0[indices] <- fit$O0
        O1[indices] <- fit$O1
        score <- score + fit$score
        information <- information + fit$information
        bridge_error <- bridge_error + fit$bridge_error
    }
    list(O0 = O0, O1 = O1, score = score, information = information, bridge_error = bridge_error)
}

