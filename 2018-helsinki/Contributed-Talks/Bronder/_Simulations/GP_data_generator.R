generate_GP_data <- function(N) {
    x <- runif(N, -10, 10)
    y <- x + 4 * x ^ 2 - x ^ 3 + 100 * sin(x * 2)
    y <- as.vector(((y - 133.4722) / 400.7542) + rnorm(N, 0, .1))
    N_predict <- N
    x_predict <- runif(N_predict, -10, 10)
    stan_data <- list(N = N, x = x, y = y, N_predict = N_predict, x_predict = x_predict)
}