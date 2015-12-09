library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lfe)
library(JJHmisc)

CreateClusteredData <- function(num.groups, sample.size, beta, randomize.by.group = FALSE){
    individual <- 1:sample.size # index for individuals
    group <- sample(1:num.groups, size = sample.size, replace = TRUE) # randomly assign everyone to a group
    group.effect <- rnorm(num.groups) # create a group-specific effect
    eta <- group.effect[group] # give each individual their group-specific effect
    epsilon <- rnorm(sample.size) # given each individual an individual-specific effect
    if (randomize.by.group){
        group.trt.assignment <- rbinom(num.groups, 1, 0.5) 
        trt <- group.trt.assignment[group]
    } else {
        trt <- rbinom(sample.size, 1, 0.5)
    }
    y <- beta * trt + epsilon + eta 
    data.frame(y, individual, group, trt, epsilon, eta)
}

CreateClusteredData(2, 10, 1, randomize.by.group = FALSE)


SimExperiment <- function(num.groups, sample.size, beta, randomize.by.group = FALSE){
    "Simulate running an analysis of the experiment with linear regression"
    df <- CreateClusteredData(num.groups, sample.size, beta, randomize.by.group = randomize.by.group)
    m <- lm(y ~ trt, data = df)
    c(as.numeric(coef(m)[2]), as.numeric(sqrt(diag(vcov(m))[2])))
}

num.replications <- 1000
df.sim.results <- replicate(num.replications, SimExperiment(20, 1000, 1)) %>%
    t %>%
    as.data.frame %>%
    set_colnames(c("beta.hat", "se"))


# median value of standard errors
(median.se <- df.sim.results %$% se %>% median)

# standard deviation of the estimated betas
(se <- df.sim.results %$% beta.hat %>% sd)

# CI comparison
alpha <- 0.05
(empirical.ci <- df.sim.results %$% beta.hat %>% quantile(., probs = c(alpha/2, 1-alpha/2)))
(calculated.ci <- c(1 - 1.96 * median.se, 1 + 1.96 * median.se))

g <- ggplot(data = df.sim.results, aes(x = beta.hat)) +
    geom_density() +
    theme_bw() +
    xlab("Simulated estimates of beta") +
    geom_vline(xintercept = empirical.ci, colour = "red") +
    geom_vline(xintercept = calculated.ci, colour = "blue")

JJHmisc::writeImage(g, "comparison", path = "", width = 5, height = 5)

############################################
## What about not accounting for clustering? 
############################################

num.replications <- 1000
df.sim.results <- replicate(num.replications, SimExperiment(20, 1000, 1, TRUE)) %>%
    t %>%
    as.data.frame %>%
    set_colnames(c("beta.hat", "se"))

(median.se <- df.sim.results %$% se %>% median)
(se <- df.sim.results %$% beta.hat %>% sd)


# What fraction of beta.hat are in a 95\% CI? 
alpha <- 0.05
empirical.ci <- df.sim.results %$% beta.hat %>% quantile(., probs = c(alpha/2, 1-alpha/2))
calculated.ci <- c(1 - 1.96 * median.se, 1 + 1.96 * median.se)

g <- ggplot(data = df.sim.results, aes(x = beta.hat)) +
    geom_density() +
    theme_bw() +
    xlab("Simulated estimates of beta") +
    geom_vline(xintercept = empirical.ci, colour = "red") +
    geom_vline(xintercept = calculated.ci, colour = "blue")

JJHmisc::writeImage(g, "no_cluster", path = "", width = 5, height = 5)

######################
# Random Effects model
######################


SimExperimentRE <- function(num.groups, sample.size, beta, randomize.by.group = FALSE){
    "Simulate running an analysis of the experiment with linear regression"
    df <- CreateClusteredData(num.groups, sample.size, beta, randomize.by.group = randomize.by.group)
    m <- lmer(y ~ trt + (1|group), data = df)
    c(as.numeric(fixef(m)[2]), as.numeric(sqrt(diag(vcov(m))[2])))
}

num.replications <- 1000
df.sim.results <- replicate(num.replications, SimExperimentRE(20, 1000, 1, TRUE)) %>%
    t %>%
    as.data.frame %>%
    set_colnames(c("beta.hat", "se"))

# what's the median standard error value? 
(median.se <- df.sim.results %$% se %>% median)

# what's the standard deviation? 
(se <- df.sim.results %$% beta.hat %>% sd)

# What fraction of beta.hat are in a 95\% CI? 
alpha <- 0.05
empirical.ci <- df.sim.results %$% beta.hat %>% quantile(., probs = c(alpha/2, 1-alpha/2))
calculated.ci <- c(1 - 1.96 * median.se, 1 + 1.96 * median.se)

g <- ggplot(data = df.sim.results, aes(x = beta.hat)) +
    geom_density() +
    theme_bw() +
    xlab("Simulated estimates of beta") +
    geom_vline(xintercept = empirical.ci, colour = "red") +
    geom_vline(xintercept = calculated.ci, colour = "blue")

JJHmisc::writeImage(g, "random_effects", path = "", width = 5, height = 5)

######################
# Fixed Effects Model 
######################


SimExperimentFE <- function(num.groups, sample.size, beta, randomize.by.group = FALSE){
    "Simulate running an analysis of the experiment with linear regression"
    df <- CreateClusteredData(num.groups, sample.size, beta, randomize.by.group = randomize.by.group)
    m <- felm(y ~ trt |0|0|group, data = df)
    c(as.numeric(coef(m)[2]), as.numeric(m$cse[2]))
}


num.replications <- 1000
df.sim.results <- replicate(num.replications, SimExperimentRE(20, 1000, 1, TRUE)) %>%
    t %>%
    as.data.frame %>%
    set_colnames(c("beta.hat", "se"))

# what's the median standard error value? 
(median.se <- df.sim.results %$% se %>% median)
# what's the standard deviation? 
(se <- df.sim.results %$% beta.hat %>% sd)

# What fraction of beta.hat are in a 95\% CI? 
alpha <- 0.05
empirical.ci <- df.sim.results %$% beta.hat %>% quantile(., probs = c(alpha/2, 1-alpha/2))
calculated.ci <- c(1 - 1.96 * median.se, 1 + 1.96 * median.se)

g <- ggplot(data = df.sim.results, aes(x = beta.hat)) +
    geom_density() +
    theme_bw() +
    xlab("Simulated estimates of beta") +
    geom_vline(xintercept = empirical.ci, colour = "red") +
    geom_vline(xintercept = calculated.ci, colour = "blue")

JJHmisc::writeImage(g, "clustering_fix", path = "", width = 5, height = 5)

print(g)
