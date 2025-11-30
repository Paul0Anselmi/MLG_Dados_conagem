##########################################################
# Trabalho: DADOS MDVIS - Número de visitas ao médico

##########################################################
### LIBS
library(doBy)
library(MASS)

# Leitura dos dados e infos iniciais
mdvis <- read.csv(
  'https://vincentarelbundock.github.io/Rdatasets/csv/COUNT/mdvis.csv'
)

str(mdvis)
names(mdvis)
summary(mdvis)

# remover a coluna rownames
mdvis$rownames <- NULL   
str(mdvis)
names(mdvis) #### confirma que foi removida a coluna, ufa


# ANÁLISE DESCRITIVA
##########################################################

# Histograma da resposta (contagem)
hist(mdvis$numvisit,
     main = 'Número de visitas ao médico (3 meses)',
     xlab = 'numvisit')

# Tabela de frequências
table(mdvis$numvisit)

# Média e variância da contagem
mean(mdvis$numvisit)
var(mdvis$numvisit)





tab_reform  <- summaryBy(numvisit ~ factor(reform),
                         data = mdvis, FUN = c(mean, sd))
tab_agegrp  <- summaryBy(numvisit ~ factor(agegrp),
                         data = mdvis, FUN = c(mean, sd))
tab_educ    <- summaryBy(numvisit ~ factor(educ),
                         data = mdvis, FUN = c(mean, sd))
tab_badh    <- summaryBy(numvisit ~ factor(badh),
                         data = mdvis, FUN = c(mean, sd))

# 2) Padronizar nomes e adicionar a qual variável pertence 

names(tab_reform)[1] <- 'categoria'
tab_reform$variavel  <- 'reform'

names(tab_agegrp)[1] <- 'categoria'
tab_agegrp$variavel  <- 'agegrp'

names(tab_educ)[1]   <- 'categoria'
tab_educ$variavel    <- 'educ'

names(tab_badh)[1]   <- 'categoria'
tab_badh$variavel    <- 'badh'

# 3) Empilhar tudo em uma tabela só 

tabela_resumo <- rbind(tab_reform, tab_agegrp, tab_educ, tab_badh)

tabela_resumo <- tabela_resumo[, c('variavel', 'categoria',
                                   'numvisit.mean', 'numvisit.sd')]

tabela_resumo


# AJUSTE POISSON (LINK LOG)
##########################################################

aj1_mdvis <- glm(
  numvisit ~ reform + factor(educ) + factor(agegrp),
  family = poisson(link = 'log'),
  data   = mdvis
)

summary(aj1_mdvis)

# Razões de taxa (exp(beta))
exp(coefficients(aj1_mdvis))

# VERIFICAÇÃO DE SOBREDISPERSÃO
##########################################################

# Resíduos de Pearson
rp <- residuals(aj1_mdvis, type = 'pearson')

# Estimativa de phi (fator de dispersão)
phi_hat <- sum(rp^2) / aj1_mdvis$df.residual
phi_hat

###### não se ajustou bem a Poisson, por isso vamos seguir para a binomial negativa. o Valor do phi_hat deu 5.749279

# AJUSTE BINOMIAL NEGATIVA (glm.nb)
##########################################################

aj_cheio <- glm.nb(
  numvisit ~ reform + badh + loginc +
    factor(educ) + factor(agegrp),
  link = log,
  data = mdvis
)

summary(aj_cheio)
exp(coef(aj_cheio))
AIC(aj_cheio)

################ VAMOS SEGUIR COM A BIN NEG


aj_step <- stepAIC(
  aj_cheio,
  direction = 'both',   # ou 'backward'
  trace = TRUE          # mostra o caminho
)

summary(aj_step)
exp(coef(aj_step))      # razões de taxa do modelo final
AIC(aj_cheio, aj_step)

aj_final <- aj_step
summary(aj_final)
exp(coef(aj_final))

### dispersão no modelo NB
rp_final <- resid(aj_final, type = 'pearson')
phi_nb   <- sum(rp_final^2) / aj_final$df.residual
phi_nb

### Matriz X, pesos, alavanca (h)

X <- model.matrix(aj_final)
w <- aj_final$weights
W <- diag(w)

auxh <- solve(t(X) %*% W %*% X)
H    <- sqrt(W) %*% X %*% auxh %*% t(X) %*% sqrt(W)
h    <- diag(H)

plot(fitted(aj_final), h,
     xlab = 'Valores ajustados',
     ylab = 'Alavanca')

identify(fitted(aj_final), h, n = 3)

### Resíduos deviance ajustados

aux_tdi <- resid(aj_final, type = 'deviance')
tdi     <- aux_tdi * sqrt(1 / (1 - h))

plot(fitted(aj_final), tdi,
     xlab = 'Valores ajustados',
     ylab = 'Resíduo deviance ajustado')
abline(h = c(-3, -2, 0, 2, 3),
       lty = c(2, 3, 1, 3, 2), col = 'grey')


# Residuos de person ajustados

aux_tsi <- resid(aj_final, type = 'pearson')
tsi     <- aux_tsi * sqrt(1 / (1 - h))

plot(fitted(aj_final), tsi,
     xlab = 'Valores ajustados',
     ylab = 'Resíduo de Pearson ajustado')
abline(h = c(-3, -2, 0, 2, 3),
       lty = c(2, 3, 1, 3, 2), col = 'grey')


# comparar tdi x tsi
plot(tdi, tsi,
     xlab = 'Resíduo deviance ajustado',
     ylab = 'Resíduo de Pearson ajustado')
abline(0, 1, lty = 2)



#### Distancia de cook

ldi <- h * (tsi^2) / (1 - h)

plot(ldi,
     xlab = 'Índice',
     ylab = 'Distância de Cook',
     type = 'h')


fit.model <- aj_final

source('C:/Users/Arthur Scheffer/Downloads/envelope_negbin.R')
