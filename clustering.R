library(magrittr)
library(dplyr)
library(ggplot2)


SimExperiment <- function(num.groups, sample.size, beta){
    individuals <- 1:sample.size
    group <- sample(1:num.groups, size = sample.size, replace = TRUE)
    group.effects <- rnorm(num.groups)
    individual.effects <- rnorm(sample.size)
    trt.assignment <- rbinom(sample.size, 1, 0.5)
    df.individual <- data.frame(individual.effect = individual.effects,
                                group = group,
                                trt = trt.assignment)
    df.group.effects <- data.frame(group = 1:num.groups, group.effect = group.effects)
    ## simulate experiment
    df.final <- df.individual %>%
    left_join(df.group.effects, by = "group") %>%
        mutate(y = beta * trt + individual.effect + group.effect)
    m <- lm(y ~ trt, data = df.final)
    c(as.numeric(coef(m)[2]), as.numeric(sqrt(diag(vcov(m)))[2]))
}

num.groups <- 10
sample.size <- 1000
beta <- 1

SimExperiment(10, 1000, 1)

sim.results <- lapply(1:1000, function(x) SimExperiment(10, 1000, 1))

df <- as.data.frame(do.call("rbind", sim.results))
colnames(df) <- c("effect", "se")

ggplot(data = df, aes(x = se)) + geom_density()

## 
mean(df$se)

mean(df$effect)

sd(df$effect)
