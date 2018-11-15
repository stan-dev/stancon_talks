library(deSolve)
library(tidyverse)

coag.ode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    cascade_delay <- 1/(1+exp(c-b*t))
    tfpi <- (1-tfpi_min)/(1+exp(-ct+bt*t)) + tfpi_min

    r_FIIa <- k_FIIa*cascade_delay*tfpi*((FII/g_FII)^c_FIIa)/(K_FIIa + (FII/g_FII)^c_FIIa)
    r_AT <- k_AT*FIIa*(AT/g_AT)
    r_clot <- k_clot*FIIa*((Fg/9e-6)^c_clot)/(K_clot + (Fg/9e-6)^c_clot)
    r_lys <- k_lys*tPA*(Fn^c_lys)/(K_lys + Fn^c_lys)
    r_PAI <- k_PAI*tPA*PAI
    
    dFII <- g_FII*(-r_FIIa)#-r_prothrombinase)                 # pct activity
    dFIIa <- r_FIIa -r_AT           # 10s of nmol/L
    dAT <- -g_AT*r_AT            # pct activity
    dFg <- -r_clot                  # 1 mg/dl = 29.41 nmol/L
    dFn <- 1e7*(r_clot - r_lys)           # mm of clot
    dtPA <- -r_PAI                  # 1 ng/mL = 14.29 pmol/L     
    dPAI <- -r_PAI                  # 1 ng/mL = 23.26 pmol/L
    
    list(c(dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI))
  })
}

parameters <- c(g_FII = 100/1.4e-6, g_AT = 100/3.4e-6,
                c = 5.0, b = 0.4e-1,
                ct = 10, bt = 0.15e-1, tfpi_min = 0.2,
                k_FIIa = 3.5e-9, c_FIIa = 1, K_FIIa = 1.4e-6,
                #k_prothrombinase = 1e8, prothrombinase_max = 1.5e-10, ctt = 7, btt = 6e-3,
                k_AT = 1.6e4,
                k_clot = 3.0, c_clot = 1, K_clot = 0.75,
                k_lys = 1, c_lys = 1, K_lys = 0.5,
                k_PAI = 4.5e5)
state      <- c(FII = 1e2, FIIa = 0, AT = 100, Fg = 9e-6, Fn = 0, tPA = 7e-11, PAI = 4e-10)
times      <- seq(0, 1800, by = 30)

out <- ode(y = state, times = times, func = eight, parms = parameters, method = "bdf",  atol = 1e-6, rtol = 1e-5) %>% as.data.frame %>% as_tibble

out %>%
  gather(state, value, -time) %>%
  mutate(state = factor(state, levels = c("FII", "FIIa", "AT", "Fg", "Fn", "tPA", "PAI"))) %>%
  ggplot(aes(time,value)) + geom_line() + facet_grid(state ~ ., scales = "free") +
  geom_line(aes(time,value), color = "red", data = mit2)

