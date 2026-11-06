# Run simulations or replay stored simulation records.
options(stringsAsFactors = FALSE)
invisible(try(suppressWarnings(Sys.setlocale('LC_ALL', 'Chinese (Simplified)_China.utf8')), silent = TRUE))
Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
root <- normalizePath(file.path(dirname(script), '..'), winslash = '/')
for (file in c('data.R', 'methods.R', 'chiang_loo.R', 'experiment.R')) source(file.path(root, 'code', file))
arg <- function(k, d) {
    value <- grep(paste0('^--', k, '='), commandArgs(TRUE), value = TRUE)
    if (length(value)) sub('^[^=]*=', '', value[1]) else d
}
mode <- arg('mode', 'replay'); job <- as.integer(arg('job', '0'))
groups <- read.csv(file.path(root, 'code/settings.csv'))
if (job) groups <- groups[groups$job == job, ]
stopifnot(nrow(groups) > 0L, mode %in% c('replay', 'run'))
out <- arg('output', file.path(root, 'results')); dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (mode == 'replay') {
    checks <- list()
    for (i in seq_len(nrow(groups))) {
        g <- groups[i, ]; saved <- readRDS(file.path(root, 'results/raw', paste0(g$id, '.rds')))
        B <- dim(saved$raw)[1]
        ids <- unique(c(1L, 37L, B, head(which(apply(!is.finite(saved$parts[, , 'v', drop = FALSE]), 1, any)), 1L)))
        for (b in ids) {
            actual <- compute_trial(b + g$offset, g)
            for (method in c('lr', 'ols', 'dof', 'loo', 'cloo', 'ploo', 'oracle')) {
                fields <- paste0(c('u_', 'v_'), method)
                a <- as.numeric(actual$raw[, fields]); s <- as.numeric(saved$raw[b, , fields])
                stopifnot(identical(is.finite(a), is.finite(s)))
                k <- is.finite(s); error <- if (any(k)) max(abs(a[k] - s[k])/(1 + abs(s[k]))) else 0
                stopifnot(error < 1e-8)
                checks[[length(checks) + 1L]] <- data.frame(group = g$id, trial = b + g$offset,
                    method = method, max_relative_error = error, finite_pattern_matches = TRUE)
            }
        }
        cat('REPLAYED', g$id, '\n')
    }
    write.csv(do.call(rbind, checks), file.path(out, if (job) paste0('code_check_job', job, '.csv') else 'code_check.csv'), row.names = FALSE)
} else {
    B <- as.integer(arg('B', '10000')); stopifnot(B > 0L)
    dir.create(file.path(out, 'raw'), showWarnings = FALSE)
    for (i in seq_len(nrow(groups))) {
        g <- groups[i, ]; dest <- file.path(out, 'raw', paste0(g$id, '.rds'))
        if (file.exists(dest)) { cat('EXISTS', g$id, '\n'); next }
        raw <- parts <- pooled_components <- meta <- NULL
        for (b in seq_len(B)) {
            fit <- compute_trial(b + g$offset, g)
            if (is.null(raw)) {
                raw <- array(NA_real_, c(B, dim(fit$raw)), dimnames = c(list(as.character(seq_len(B))), dimnames(fit$raw)))
                parts <- array(NA_real_, c(B, dim(fit$parts)), dimnames = c(list(as.character(seq_len(B))), dimnames(fit$parts)))
                pooled_components <- array(NA_real_, c(B, dim(fit$pooled_components)),
                    dimnames = c(list(as.character(seq_len(B))), dimnames(fit$pooled_components)))
                meta <- matrix(NA_real_, B, length(fit$meta), dimnames = list(NULL, names(fit$meta)))
            }
            raw[b, , ] <- fit$raw; parts[b, , ] <- fit$parts; pooled_components[b, , , ] <- fit$pooled_components; meta[b, ] <- fit$meta
            if (b %% 1000L == 0L) { cat(g$id, b, '/', B, '\n'); flush.console() }
        }
        connection <- xzfile(dest, open = 'wb', compression = 9)
        tryCatch(saveRDS(list(group = g, raw = raw, parts = parts, pooled_components = pooled_components,
            meta = meta, offset = g$offset), connection), finally = close(connection))
    }
}
