# Reproduce the ACTG175 Z3 application with three baseline adjustment sets.
options(stringsAsFactors = FALSE, digits = 16)
invisible(try(suppressWarnings(Sys.setlocale('LC_ALL', 'Chinese (Simplified)_China.utf8')), silent = TRUE))
suppressPackageStartupMessages(library(survival))
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
root <- normalizePath(file.path(dirname(script), '..', '..'), winslash = '/')
arg <- function(name, default) {
    value <- grep(paste0('^--', name, '='), commandArgs(TRUE), value = TRUE)
    if (length(value)) sub('^[^=]*=', '', value[1]) else default
}
out <- arg('output', file.path(root, 'application/results'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (file in c('data.R', 'methods.R', 'chiang_loo.R')) source(file.path(root, 'code', file))

raw <- read.table(file.path(root, 'application/data/ACTG175.txt'), header = TRUE)
selected <- raw$arms %in% c(0, 3) & raw$strat == 3
d <- raw[selected, , drop = FALSE]
A <- as.integer(d$arms == 3)
n <- nrow(d)
q <- mean(A)
sets <- list(
    p5 = c('cd40', 'preanti', 'cd80', 'age', 'wtkg'),
    p8 = c('cd40', 'preanti', 'cd80', 'age', 'wtkg', 'karnof', 'symptom', 'gender'),
    p14 = c('cd40', 'preanti', 'cd80', 'age', 'wtkg', 'karnof', 'symptom', 'gender',
            'hemo', 'homo', 'drugs', 'race', 'oprior', 'z30')
)
stopifnot(n == 434L, sum(A) == 221L, sum(d$cens) == 148L,
          sum(d$cens[A == 1]) == 56L, sum(d$cens[A == 0]) == 92L,
          !anyNA(d[, c(sets$p14, 'days', 'cens')]))

der <- compute_derived_outcomes(d$days, d$cens, A)
lr <- fit_logrank(d$days, d$cens, A)
checks <- list(data.frame(set = 'all', check = 'derived_score_bridge', error = abs(der$bridge_error)))
test_rows <- diagnostics <- list()
for (sn in names(sets)) {
    X <- scale(as.matrix(d[sets[[sn]]]), TRUE, TRUE)
    p <- ncol(X)
    inv <- mf_inverse(X)
    H <- tcrossprod(X %*% inv$W, X)
    h <- diag(H)
    arm_rank <- vapply(0:1, function(a) qr(cbind(1, X[A == a, , drop = FALSE]))$rank, integer(1))
    stopifnot(inv$rank == p, all(arm_rank == p + 1L))
    pooled <- pooled_fit(X, A, der$O0, der$O1, der$information/n, methods = c('cloo', 'ploo'))
    robin <- shared_robin(X, A, der$O0, der$O1, der)
    chiang <- cmo_grid(cbind(1, X), A, der$O0, der$O1, p + 1)$parts[1, ]
    fits <- rbind(
        lr = lr[c('u', 'v')],
        ols = c(u = unname(robin['u_ols']), v = unname(robin['v_ols'])),
        dof = c(u = unname(robin['u_dof']), v = unname(robin['v_dof'])),
        loo = c(u = unname(chiang['u_check']), v = unname(chiang['mHC3'])),
        pooled[, c('u', 'v')]
    )
    for (method in rownames(fits)) {
        u <- fits[method, 'u']
        v <- fits[method, 'v']
        test_rows[[length(test_rows) + 1L]] <- data.frame(
            set = sn, p, n, events = sum(d$cens), method,
            u, sqrt_n_u = sqrt(n) * u, variance = v, sigma = sqrt(v),
            se_u = sqrt(v/n), statistic = sqrt(n) * u/sqrt(v),
            I = if (method %in% rownames(pooled)) pooled[method, 'I'] else NA_real_,
            II = if (method %in% rownames(pooled)) pooled[method, 'II'] else NA_real_
        )
    }
    diagnostics[[sn]] <- data.frame(set = sn, p, n, rank0 = arm_rank[1], rank1 = arm_rank[2],
                                   max_h = max(h), fields = paste(sets[[sn]], collapse = ';'))

    # Dense off-diagonal formulas independently check the pooled computations.
    y <- ifelse(A == 1, der$O1/q, der$O0/(1 - q))
    response <- der$O0 + der$O1
    off <- row(H) != col(H)
    u_literal <- mean(A * der$O1 - (1 - A) * der$O0 + (q - A) * (drop(H %*% y) - h * y))
    I_literal <- der$information/n - q * (1 - q)/n *
        sum(((1 + outer(h, h, '+') + H) * H * tcrossprod(response))[off])
    K1 <- H * outer(der$O1, der$O1, '+')
    K0 <- H * outer(der$O0, der$O0, '+')
    II_literal <- sum(((1 - q) * K1^2 + q * K0^2)[upper.tri(H)])/n
    checks[[length(checks) + 1L]] <- data.frame(set = sn, check = 'ploo_dense',
        error = max(abs(c(u_literal - pooled['ploo', 'u'],
                          I_literal - pooled['ploo', 'I'], II_literal - pooled['ploo', 'II']))))
    C <- sweep(H, 1, 1 - h, '/')
    diag(C) <- 0
    u_literal <- mean(A * der$O1 - (1 - A) * der$O0 + (q - A) * drop(C %*% y))
    I_literal <- der$information/n - q * (1 - q)/n *
        sum(((2 * C - crossprod(C) + C * t(C)) * tcrossprod(response))[off])
    K1 <- sweep(C, 2, der$O1, '*')
    K0 <- sweep(C, 2, der$O0, '*')
    II_literal <- sum(((1 - q) * (K1 + t(K1))^2 + q * (K0 + t(K0))^2)[upper.tri(H)])/n
    checks[[length(checks) + 1L]] <- data.frame(set = sn, check = 'cloo_dense',
        error = max(abs(c(u_literal - pooled['cloo', 'u'],
                          I_literal - pooled['cloo', 'I'], II_literal - pooled['cloo', 'II']))))
}
tests <- do.call(rbind, test_rows)
checks <- do.call(rbind, checks)
stopifnot(nrow(tests) == 18L, all(tests$variance > 0), all(is.finite(tests$statistic)),
          max(checks$error) < 1e-10)
write.csv(tests, file.path(out, 'test_statistics.csv'), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, 'design_diagnostics.csv'), row.names = FALSE)
write.csv(checks, file.path(out, 'formula_checks.csv'), row.names = FALSE)
write.csv(data.frame(source_row = which(selected), pidnum = d$pidnum),
          file.path(out, 'sample_index.csv'), row.names = FALSE)
write.csv(data.frame(arm = c('zidovudine', 'didanosine'), n = c(sum(A == 0), sum(A == 1)),
                    events = c(sum(d$cens[A == 0]), sum(d$cens[A == 1]))),
          file.path(out, 'sample_summary.csv'), row.names = FALSE)

display <- data.frame(p = tests$p, method = tests$method,
                      numerator = sprintf('%.4f', tests$sqrt_n_u),
                      sigma = sprintf('%.4f', tests$sigma), statistic = sprintf('%.3f', tests$statistic))
write.csv(display, file.path(out, 'table.csv'), row.names = FALSE)
table <- c('\\begin{table}[htbp]', '\\centering\\small',
           '\\caption{Results for reanalyzing the ACTG175 subgroup.}',
           '\\label{tab:actg175-application}', '\\setlength{\\tabcolsep}{5pt}',
           '\\renewcommand{\\arraystretch}{1.05}', '\\begin{tabular}{@{}rrrrrrr@{}}', '\\toprule')
fields <- c('numerator', 'sigma', 'statistic')
titles <- c('Numerator $\\sqrt{n}\\widehat U$',
            'Estimated standard deviation $\\widehat\\sigma$',
            'Test statistic $\\sqrt{n}\\widehat U/\\widehat\\sigma$')
for (j in seq_along(fields)) {
    if (j > 1) table <- c(table, '\\midrule')
    header <- c('$p$', paste0('$\\widehat\\chi_{\\mathrm{', c('lr', 'ols', 'dof', 'loo', 'cloo', 'ploo'), '}}$'))
    table <- c(table, paste0('\\multicolumn{7}{c}{', titles[j], '} \\\\'), '\\midrule',
               paste0(paste(header, collapse = ' & '), ' \\\\'), '\\midrule')
    for (p in c(5, 8, 14)) {
        values <- display[display$p == p, fields[j]]
        table <- c(table, paste0(paste(c(p, values), collapse = ' & '), ' \\\\'))
    }
}
writeLines(c(table, '\\bottomrule', '\\end{tabular}', '\\end{table}'), file.path(out, 'table.tex'))
cat('Saved ACTG175 Z3 results: 434 patients, 148 events, 18 method settings.\n')
print(tests[, c('p', 'method', 'sqrt_n_u', 'sigma', 'statistic')], row.names = FALSE)
