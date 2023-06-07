pacman::p_load(readr, tidyverse, refund,lubridate, refund.shiny)

data("DTI")

DTI2 <- DTI %>% 
  drop_na() %>% 
  filter(visit == 1 & case == 1)

#

passat <- DTI2$pasat

tract <- 1:93

cca <- DTI2$cca

#

fpca_res <- fpca.sc(cca, argvals = tract, pve = 0.97)

plot_shiny(fpca_res)

m <- length(tract)
efn <- fpca_res$efunctions*sqrt(m)
eval <- fpca_res$evalues/m
scr <- fpca_res$scores/sqrt(m)
npc <- fpca_res$npc

cbind(tract,efn) %>% 
  as_tibble() %>% 
  rename(
    pc1 = V2,
    pc2 = V3,
    pc3 = V4,
    pc4 = V5,
    pc5 = V6,
    pc6 = V7,
    pc7 = V8,
    ) %>% 
  pivot_longer(cols = starts_with("pc")) %>% 
  group_by(tract) %>% 
  ggplot(aes(tract,value,group = name, color = name))+
  geom_line()

#

k.pc <- 1
effect <- sqrt(eval[k.pc])*efn[,k.pc]
mu_hat <- fpca_res$mu


cbind(tract,efn[,k.pc]) %>% 
  as_tibble() %>% 
  ggplot(aes(tract, V2))+
  geom_line()+
  ylim(c(-2,2))


mu_hat %>% 
  as_tibble_col(column_name = "mu.hat") %>% 
  mutate(
    effect_Plus = mu.hat + effect,
    effect_less = mu.hat - effect,
    tract = tract
  ) %>% 
  pivot_longer(cols = 1:3) %>% 
  group_by(tract) %>% 
  ggplot(aes(tract, value, group = name, color = name))+
  geom_line()

#
out = lm(passat ~ scr) ## Multiple linear regression
# summary(out)
beta_hat = out$coefficients
beta_hat

summary(out)

beta_fn_hat = efn %*% as.matrix(beta_hat[-1], col = 1) 

cbind(tract,beta_fn_hat) %>% 
  as_tibble() %>% 
  ggplot(aes(tract,V2))+
  geom_line()+
  ylim(c(-2000,800))



set.seed(12)
n.crv <- 3
n <- nrow(cca)
sel.crv <- sample(1:n, size=n.crv, replace = FALSE)


yhat <- t(fpca_res$Yhat[sel.crv,])


cbind(tract,yhat) %>% 
  as_tibble() %>% 
  pivot_longer(cols = contains("_")) %>% 
  group_by(tract) %>% 
  ggplot(aes(tract,value,group = name,color=name))+
  geom_line()



par(mfrow=c(3,3))
for(i in 1:3){
  ind <- sel.crv[i]
  demeaned <- fpca_res$Yhat[ind,]-as.vector(fpca_res$mu)
  
  matplot(tract, t(fpca_res$Yhat[sel.crv,]-t(matrix(rep(fpca_res$mu,3), nrow=93))), 
          type='l', lwd=2, lty=1, col = 'light grey',
          xlab="miles", ylab="speed (demeaned)", main="")
  lines(tract, demeaned, type='l', lwd=2, col='red')
  
  
  plot(tract, beta_fn_hat, type='l', lwd=2,
       xlab="miles", ylab = "estimated coefficient fn", main="")
  plot(tract, demeaned*beta_fn_hat,type='l', lwd=2, col='blue',
       xlab="miles", ylab = "", ylim=c(-55, 70),
       main=round(mean(demeaned*beta_fn_hat), 2))
}


par(mfrow=c(1,1))
plot(passat, out$fitted, cex=0.5, ylab="Fitted", xlab="Observed")
abline(a = 0, b = 1)


Rsq = 1-sum((out$residuals)^2)/sum((passat- mean(passat))^2)
Rsq



####
X <- as.matrix(cca) # functional covariate
Y <- passat    # scalar response


myDat <- data.frame(X, Y)

fit <- pfr(Y ~ lf(X, k = 10, bs = "cr"), method = "REML", data = myDat)
coef <- coef(fit)


cbind(coef$X.argvals, coef$value) %>% 
  as_tibble() %>% 
  ggplot(aes(V1,V2))+
  geom_point()


# -----------------------------------

par(mfrow=c(3,3))

beta_fn_hat0 <- beta_fn_hat  # saving beta(t) from fPCA approach

fit <- pfr(Y ~ lf(X, k = 10, bs = "cr"), method = "REML", data = myDat)
coef <- coef(fit)
beta_fn_hat <- coef$value

for(i in 1:3){
  ind <- sel.crv[i]
  demeaned <- fpca_res$Yhat[ind,]-as.vector(fpca_res$mu)
  
  matplot(tract, t(fpca_res$Yhat[sel.crv,]-t(matrix(rep(fpca_res$mu,3), nrow=93))), 
          type='l', lwd=2, lty=1, col = 'light grey',
          xlab="miles", ylab="speed (demeaned)", main="")
  lines(tract, demeaned, type='l', lwd=2, col='red')
  
  
  plot(tract, beta_fn_hat, type='l', lwd=2,
       xlab="miles", ylab = "estimated coefficient fn", main="")
  plot(tract, demeaned*beta_fn_hat,type='l', lwd=2, col='blue',
       xlab="miles", ylab = "", main=round(mean(demeaned*beta_fn_hat), 2)) 
  
}


par(mfrow=c(1,1))
plot(rowSums( sapply(1:7, function(a) (t(efn) %*% beta_fn_hat/93)[a] * efn[,a]) ) , type='l', lwd=2, col = "red", ylab="", xlab="miles")
lines(tract, beta_fn_hat0, type='l', lwd=2)





#
fpca_res <- fpca.sc(cca, argvals = tract, pve = 0.97)
Xhat <- fpca_res$Yhat
Yhat <- predict(fit, newdata = list(X = Xhat))

# goodness-of-fit
par(mfrow=c(1,1))
plot(Y, Yhat, cex=0.5, ylab="Fitted", xlab="Observed")
abline(a = 0, b = 1)

Rsq = 1-sum((Y- as.vector(Yhat))^2)/sum((Y - mean(Y))^2)
Rsq

####################################################################
###################################################################
library(tidymodels)

marathon_results_2017 <-readr::read_csv("marathon_results_2017.csv")

glimpse(marathon_results_2017)

quartis <- quantile(marathon_results_2017$Age, probs = c(0, 0.25, 0.5, 0.75, 1))


marathon_df<- marathon_results_2017 %>% 
  dplyr::select(-c(`...1`,Bib,City,State,Country,
                   Citizen,`...10`,Division,Gender,
                   Overall,`Proj Time`)) %>% 
  janitor::clean_names() %>% 
  mutate(
    age_group = cut(age,
                    breaks = quartis, labels = c("Q1", "Q2", "Q3", "Q4"),
                    include.lowest = TRUE)
         ) %>% 
  rename(gender = m_f)


# set.seed(369)
# marathon_df2 <- marathon_df %>%
#   group_by(gender, age_group) %>%
#   slice_sample(n = 25, replace = FALSE) %>% 
#   ungroup()


marathon_df2<- marathon_df %>% 
  filter(gender == "M") %>% 
  group_by(age_group) %>%
  slice_sample(n = 50, replace = FALSE) %>% 
  ungroup()



mutate_seconds <- function(x){
  lubridate::hour(x)* 3600 + lubridate::minute(x) * 60 + lubridate::second(x) *1
}
########################################################
# Definindo o tamanho do vetor
tamanho <- 200

# Definindo os limites inferior e superior
limite_inferior <- 14.5
limite_superior <- 17.5

# Gerando o vetor de números aleatórios
#vetor_aleatorio <- runif(tamanho, limite_inferior, limite_superior)

# Arredondando os valores para inteiros
#vetor_aleatorio <- round(vetor_aleatorio)
#######################################################




marathon_df2 <- 
  marathon_df2 %>% 
  drop_na() %>% 
  mutate(
    across(starts_with(c("x","half")), lubridate::hms),
    across(starts_with(c("x","half","pace", "official_time")), mutate_seconds),
   `2` = runif(n=200, limite_inferior, limite_superior),
     `5` = (5000/ `x5k`)*3.6,
    `10` = (10000/ `x10k`)*3.6,
    `15` = (15000/ `x15k`)*3.6,
    `20` = (20000/ `x20k`)*3.6,
    `21` = (21000/ half)*3.6,
    `25` = (25000/ `x25k`)*3.6,
    `30` = (30000/ `x30k`)*3.6,
    `35` = (35000/ `x35k`)*3.6,
    `40` = (40000/ `x40k`)*3.6,
    `42` = (42000/ official_time)*3.6
  ) %>% 
  rename(
    names = name
  )

####################################
marathon_df2 %>% 
  pivot_longer(cols = 16:26) %>% 
  mutate(
    name = as.numeric(name)
  ) %>% 
  group_by(name) %>% 
  ggplot(aes(name,value, group = names, color = m_f))+
  geom_line(alpha = 0.7, color = "grey")+
  #geom_point(color = "black", alpha = 0.7)+
  ylab("Velocidade Km/h")+
  xlab("Km de maratona")+
  ggtitle("Maratona - Velocidade")
 
##

speed <- marathon_df2 %>% 
  dplyr::select(16:26) %>% 
  as.matrix()

kms <- c(2,5,10,15,20,21,25,30,35,40,42)

final_time <- marathon_df2$official_time %>% 
  as_tibble() %>% 
  as.matrix()

fpca_res <- fpca.face(Y= speed ,pve = 0.97, argvals = kms, knots = 7)

fpca_res <- fpca.sc(speed, argvals = kms, pve = 0.97, nbasis = 7)



m <- length(kms)
efn <- fpca_res$efunctions*sqrt(m)
eval <- fpca_res$evalues/m
scr <- fpca_res$scores/sqrt(m)
npc <- fpca_res$npc


cbind(kms,efn) %>% 
  as_tibble() %>% 
  rename(
    pc1 = V2,
    pc2 = V3
  ) %>% 
  pivot_longer(cols = starts_with("pc")) %>% 
  group_by(kms) %>% 
  ggplot(aes(kms,value,group = name, color = name))+
  geom_line()+
  labs(x="km",y="", title = "Autofunções estimadas" )



k.pc <- 1
effect <- sqrt(eval[k.pc])*efn[,k.pc]
mu_hat <- fpca_res$mu


cbind(kms,efn[,k.pc]) %>% 
  as_tibble() %>% 
  ggplot(aes(kms, V2))+
  geom_line()+
  ylim(c(-2,2))+
  labs(x="km",y="", title = "fCP-1" )


mu_hat %>% 
  as_tibble_col(column_name = "mu.hat") %>% 
  mutate(
    effect_Plus = mu.hat + effect,
    effect_less = mu.hat - effect,
    kms = kms
  ) %>% 
  pivot_longer(cols = 1:3) %>% 
  group_by(kms) %>% 
  ggplot(aes(kms, value, group = name, color = name))+
  geom_line()+
  labs(x="km",y="", title = "fCP-1")


out = lm(final_time ~ scr) ## Multiple linear regression
# summary(out)
beta_hat = out$coefficients
beta_hat

summary(out)


beta_fn_hat = efn %*% as.matrix(beta_hat[-1], col = 1) 

cbind(kms,beta_fn_hat) %>% 
  as_tibble() %>% 
  ggplot(aes(kms,V2))+
  geom_line()#+
  # ylim(c(-2000,800))+
  # ylab("")

set.seed(12)
n.crv <- 3
n <- nrow(speed)
sel.crv <- sample(1:n, size=n.crv, replace = FALSE)

rand_speed <- t(fpca_res$Yhat[sel.crv,])

cbind(kms,rand_speed) %>% 
  as_tibble() %>% 
  pivot_longer(cols = contains("V")) %>% 
  group_by(kms) %>% 
  ggplot(aes(kms,value,group = name,color=name))+
  geom_line()+
  ylab("")


par(mfrow=c(3,3))
for(i in 1:3){
  ind <- sel.crv[i]
  demeaned <- fpca_res$Yhat[ind,]-as.vector(fpca_res$mu)
  
  matplot(kms, t(fpca_res$Yhat[sel.crv,]-t(matrix(rep(fpca_res$mu,3), nrow=11))), 
          type='l', lwd=2, lty=1, col = 'light grey',
          xlab="miles", ylab="speed (demeaned)", main="")
  lines(kms, demeaned, type='l', lwd=2, col='red')
  
  
  plot(kms, beta_fn_hat, type='l', lwd=2,
       xlab="miles", ylab = "estimated coefficient fn", main="")
  plot(kms, demeaned*beta_fn_hat,type='l', lwd=2, col='blue',
       xlab="miles", ylab = "", ylim=c(-55, 70),
       main=round(mean(demeaned*beta_fn_hat), 2))
}
