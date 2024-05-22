Dat = read.csv("data/dataset/P2_774725.csv", header=TRUE)

M_E <- lm(Y ~ E1+E2+E3+E4, data=Dat)

M_raw <- lm(Y ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=Dat)
pdf("figures/raw_resid.pdf", width=9.5, height=5.5)
plot(resid(M_raw) ~ fitted(M_raw),
     xlab="Squared sum of IVs", ylab="Residual")
dev.off()

library(MASS)
pdf("figures/boxcox.pdf")
boxcox(M_raw, lambda=seq(-2,3,0.1))
dev.off()
# chosen by inspecting the boxcox output
p <- 2

M_trans <- lm( I(Y^p) ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=Dat)
pdf("figures/trans_resid.pdf", width=9.5, height=5.5)
plot(resid(M_trans) ~ fitted(M_trans),
     xlab="Squared sum of IVs", ylab="Residual")
dev.off()

library(leaps)
M <- regsubsets(model.matrix(M_trans)[,-1], I((Dat$Y)^p),
                nbest=1 , nvmax=5, 
                method="forward", intercept=TRUE)
temp <- summary(M)

library("knitr")
Var <- colnames(model.matrix(M_trans))
M_select <- apply(temp$which, 1, 
                  function(x) paste0(Var[x], collapse='+'))
model_summary <- kable(data.frame(cbind( model = M_select, adjR2 = temp$adjr2, BIC = temp$bic)), caption='Model Summary')

M_main <- lm( I(Y^p) ~ E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20, data=Dat)
temp <- summary(M_main)
# (Intercept)+E2+E3+E4+E3:E4+G6:G16
# E2, E3, E4, G6, and G16 chosen based on past two sections
sig_coef <- kable(temp$coefficients[ abs(temp$coefficients[,4]) <= 0.001, ], caption='Sig Coefficients')

M_2stage <- lm(I(Y^p) ~ (E2+E3+E4+G6+G16)^2, data=Dat)
temp <- summary(M_2stage)
final_coef <- kable(temp$coefficients[ abs(temp$coefficients[,3]) >= 0, ])

M_final <- lm(I(Y^p) ~ E2+E3+E4+G6:G16, data=Dat)
final_model <- summary(M_final)
