#rm(list = ls())

### Adapted from Baier etal.
# Parameters for mothers and fathers env influence is equated in order to get 
# 1 statistics out of the model

library(OpenMx)
library(here)
mxOption(NULL, "Default optimizer", "SLSQP")

estimateMCoT <- function(df, tryhard=1) {
  dat <- df
  # Fix covariates
  # Check if there are outcomes with missing covariate values. Must be removed.
  NA_sexo = is.na(dat[, c("sex_o11", "sex_o12", "sex_o21", "sex_o22")])
  NA_cohort = is.na(dat[, c("cohort_o11", "cohort_o12", "cohort_o21", "cohort_o22")])
  missy_sexo = !is.na(dat[, 5:8]) & NA_sexo
  missy_cohorto = !is.na(dat[, 5:8]) & NA_cohort
  rem_sexo = rowSums(missy_sexo) > 0
  rem_cohorto = rowSums(missy_cohorto) > 0
  # All zero so no point
  
  # Add large number for missing sex
  dat[, c("sex_o11", "sex_o12", "sex_o21", "sex_o22")][NA_sexo] = -9999
  # Add large number for missing cohort
  dat[, c("cohort_o11", "cohort_o12", "cohort_o21", "cohort_o22")][NA_cohort] = -9999
  
  
  # Filter cohorts
  dat_mx = dat
  #save(dat_mx, file=here("data-for-analysis.RData"))
  
  var_nms = colnames(dat_mx)[1:8]
  
  # Start values
  cov_ys = cov(na.omit(select(dat_mx, all_of(var_nms))))
  var_m = mean(diag(cov_ys)[1:2])
  var_f = mean(diag(cov_ys)[3:4])
  var_o = mean(diag(cov_ys)[5:8])
  cov_mf = mean(diag(cov_ys[c(3, 4), c(1, 2)]))
  reg_om = coef(lm(o11~m1+f1, dat_mx))[2]
  reg_of = coef(lm(o11~m1+f1, dat_mx))[3]
  mean_m = mean(c(dat_mx$m1, dat_mx$m2), na.rm = T)
  mean_f = mean(c(dat_mx$f1, dat_mx$f2), na.rm = T)
  
  mod_full = mxModel("mod_full",
                     
                     # Data
                     mxData(dat_mx, "raw"),
                     mxMatrix("Fu", 1, 2, F, 0, c("data.mz", "data.dz"), name = "rel"),
                     mxMatrix("Fu", 2, 1, F, 0, c("data.sex_o11", "data.sex_o12"), name = "sex_o1"),
                     mxMatrix("Fu", 2, 1, F, 0, c("data.sex_o21", "data.sex_o22"), name = "sex_o2"),
                     
                     # Parameters
                     mxMatrix("Di", 2, 2, T, sqrt(0.5 * var_m), "a1p", lbound = 0, name = "Aa1p"), # Additive genetic path parents
                     mxMatrix("Di", 2, 2, T, sqrt(0.5 * var_o), "a1o", lbound = NA, name = "Aa1o"), # Additive genetic path children shared with parents
                     mxMatrix("Di", 2, 2, T, sqrt(0.5 * var_o), "a2o", lbound = 0, name = "Aa2o"), # Additive genetic path children unique to children
                     mxMatrix("Di", 2, 2, T, sqrt(0.2 * var_m), "cp", lbound = 0, name = "Acp"), # Common environment path parents
                     mxMatrix("Di", 2, 2, T, sqrt(0.5 * var_m), "ep", lbound = 0, name = "Aep"), # Unique environment path parents
                     mxMatrix("Di", 2, 2, T, sqrt(0.4 * var_o), "eo", lbound = 0, name = "Aeo"), # Unique environment path children
                     mxMatrix("Di", 2, 2, T, sqrt(0.2 * var_o), "co", lbound = 0, name = "Aco"), # Common environment path children
                     mxMatrix("Di", 2, 2, T, reg_of, "m", name = "Abm"), # maternal effect
                     mxMatrix("Di", 2, 2, T, reg_of, "p", name = "Abf"), # Paternal effect
                     mxMatrix("Fu", 1, 1, T, cov_mf, "d", name = "Sd"), # Covariance partners (parents)
                     mxMatrix("Fu", 1, 1, T, 0.0, "w", name = "Sw"), # Correlation between A and C parents (constrained)
                     mxMatrix("Fu", 2, 1, T, mean_m, "um", name = "Um"), # Intercepts mothers
                     mxMatrix("Fu", 2, 1, T, mean_f, "uf", name = "Uf"), # Intercepts fathers
                     mxMatrix("Fu", 2, 1, T, 2, "uo", name = "Uo"), # Intercepts children
                     mxMatrix("Di", 2, 2, T, 2, "bfemaleo", name = "Bfemaleo"), # Coefficients for female children
                     mxMatrix("Fu", 1, 1, T, 2, "rho", lbound = 0, name = "rhom"), # Squared exponetial decay in common environment correlation children
                     
                     # Transformed parameters
                     mxAlgebra(a1p^2 + ep^2 + cp^2 + 2 * a1p * cp * w, name = "v2p"), # Marginal variance parents
                     mxAlgebra(a1p / sqrt(v2p), name = "za1p"), # standardized coefficients parents
                     mxAlgebra(ep / sqrt(v2p), name = "zep"),
                     mxAlgebra(cp / sqrt(v2p), name = "zcp"),
                     mxAlgebra(d / sqrt(v2p * v2p), name = "u"), # Spousal correlation
                     
                     mxAlgebra(m^2 * v2p + p^2 * v2p + 2 * m * p * d +
                                 2 * a1o * CovACo +
                                 a1o^2 + a2o^2 + co^2 + eo^2, name = "v2o"), # Marginal variance children
                     mxAlgebra(m * sqrt(v2p / v2o), name = "zm"),
                     mxAlgebra(p * sqrt(v2p / v2o), name = "zp"),
                     mxAlgebra(a1o / sqrt(v2o), name = "za1o"),
                     mxAlgebra(a2o / sqrt(v2o), name = "za2o"),
                     mxAlgebra(co / sqrt(v2o), name = "zco"),
                     mxAlgebra(eo / sqrt(v2o), name = "zeo"),
                     
                     # Genetic correlations among parent siblings
                     # aln is adjusted due to assortment, al is not.
                     mxAlgebra(rbind(1, 1/2), name = "al_pos"),
                     mxAlgebra(rbind(1, 1/2 * (1 + za1p * u * za1p +
                                                 za1p * u * zcp * w +
                                                 w * zcp * u * za1p +
                                                 w * zcp * u * zcp * w)), name = "aln_pos"),
                     mxAlgebra(rel %*% al_pos, name = "al"),
                     mxAlgebra(rel %*% aln_pos, name = "aln"),
                     
                     # Constants
                     # MOstly for manipulating matrices
                     mxMatrix("Ze", 2, 2, name = "Z"),
                     mxMatrix("Ze", 2, 1, name = "Z1"),
                     mxMatrix("Id", 2, name = "I"),
                     mxMatrix("Id", 16, name = "Ipp"),
                     mxMatrix("Di", 2, 2, F, 0.5, name = "Atr"),
                     mxMatrix("Id", 20, name = "Ioo"),
                     
                     # Constraints
                     # This is to equate the passive gene-environment correlation across generations
                     mxAlgebra(m * a1p * 0.5 +
                                 m * u * a1p * 0.5 +
                                 p * a1p * 0.5 +
                                 p * u * a1p * 0.5 +
                                 m * cp * w * 0.5 +
                                 m * u * cp * w * 0.5 +
                                 p * cp * w * 0.5 +
                                 p * u * cp * w * 0.5, name = "CovACo"), # Cov(F, A)
                     mxAlgebra(CovACo / sqrt((m^2 * v2p + p^2 * v2p + 2 * m * p * d + 1) * (1 + 1)), name = "rACo"), # Corr(F + C', A + A')
                     mxConstraint(w == rACo, name = "Equate_rae_gen"), 
                     #mxConstraint(m == p, name = "Parents_inf_equal"), 
                     
                     # Parent model
                     # -------------------------------------------------------------------------------------------
                     # S =
                     # AA |    |  
                     # -------------
                     # EA | EE |
                     # -------------
                     # CA | CE | CC
                     # There are 6 potential relationships per covariance block that may depend on parameters.
                     # wi = within person
                     # si = siblings
                     # pa = partners
                     # sp = sibling with partner (siblings in law)
                     # ps = partner with sibling (siblings in law) NB! sp & ps are not symmetric relations in the off-diagonal matrices.
                     # co = co-siblings in law
                     
                     # Here all possible correlations among parental latent variables are defined
                     # AA
                     mxAlgebra(1, name = "wiAA"),
                     mxAlgebra(aln, name = "siAA"),
                     mxAlgebra(za1p * u * za1p +
                                 za1p * u * zcp * w +
                                 w * zcp * u * za1p +
                                 w * zcp * u * zcp * w, name = "paAA"),
                     mxAlgebra(aln * za1p * u * za1p +
                                 aln * za1p * u * zcp * w +
                                 w * zcp * u * za1p +
                                 w * zcp * u * zcp * w, name = "spAA"),
                     mxAlgebra(za1p * u * za1p * aln +
                                 za1p * u * zcp * w +
                                 w * zcp * u * za1p * aln +
                                 w * zcp * u * zcp * w, name = "psAA"),
                     mxAlgebra(za1p * u * za1p * aln * za1p * u * za1p +
                                 za1p * u * za1p * aln * za1p * u * zcp * w +
                                 za1p * u * za1p * w * zcp * u * za1p +
                                 za1p * u * za1p * w * zcp * u * zcp * w +
                                 
                                 za1p * u * zcp * 1 * zcp * u * za1p +
                                 za1p * u * zcp * 1 * zcp * u * zcp * w +
                                 za1p * u * zcp * w * za1p * u * za1p +
                                 za1p * u * zcp * w * za1p * u * zcp * w +
                                 
                                 w * zcp * u * za1p * aln * za1p * u * za1p +
                                 w * zcp * u * za1p * aln * za1p * u * zcp * w +
                                 w * zcp * u * za1p * w * zcp * u * za1p +
                                 w * zcp * u * za1p * w * zcp * u * zcp * w +
                                 
                                 w * zcp * u * zcp * 1 * zcp * u * za1p +
                                 w * zcp * u * zcp * 1 * zcp * u * zcp * w +
                                 w * zcp * u * zcp * w * za1p * u * za1p +
                                 w * zcp * u * zcp * w * za1p * u * zcp * w, name = "coAA"),
                     # EE
                     mxAlgebra(1, name = "wiEE"),
                     mxAlgebra(0, name = "siEE"),
                     mxAlgebra(zep * u * zep, name = "paEE"),
                     mxAlgebra(0, name = "spEE"),
                     mxAlgebra(0, name = "psEE"),
                     mxAlgebra(zep * u * za1p * aln * za1p * u * zep +
                                 zep * u * za1p * w * zcp * u * zep +
                                 zep * u * zcp * 1 * zcp * u * zep +
                                 zep * u * zcp * w * za1p * u * zep, name = "coEE"),
                     # CC
                     mxAlgebra(1, name = "wiCC"),
                     mxAlgebra(1, name = "siCC"),
                     mxAlgebra(zcp * u * zcp +
                                 zcp * u * za1p * w +
                                 w * za1p * u * zcp +
                                 w * za1p * u * za1p * w, name = "paCC"),
                     mxAlgebra(1 * zcp * u * zcp +
                                 1 * zcp * u * za1p * w +
                                 w * za1p * u * zcp +
                                 w * za1p * u * za1p * w, name = "spCC"),
                     mxAlgebra(zcp * u * zcp * 1 +
                                 zcp * u * za1p * w +
                                 w * za1p * u * zcp * 1 +
                                 w * za1p * u * za1p * w, name = "psCC"),
                     mxAlgebra(zcp * u * zcp * 1 * zcp * u * zcp +
                                 zcp * u * zcp * 1 * zcp * u * za1p * w +
                                 zcp * u * zcp * w * za1p * u * zcp +
                                 zcp * u * zcp * w * za1p * u * za1p * w +
                                 
                                 zcp * u * za1p * aln * za1p * u * zcp +
                                 zcp * u * za1p * aln * za1p * u * za1p * w +
                                 zcp * u * za1p * w * zcp * u * zcp +
                                 zcp * u * za1p * w * zcp * u * za1p * w +
                                 
                                 w * za1p * u * zcp * 1 * zcp * u * zcp +
                                 w * za1p * u * zcp * 1 * zcp * u * za1p * w +
                                 w * za1p * u * zcp * w * za1p * u * zcp +
                                 w * za1p * u * zcp * w * za1p * u * za1p * w +
                                 
                                 w * za1p * u * za1p * aln * za1p * u * zcp +
                                 w * za1p * u * za1p * aln * za1p * u * za1p * w +
                                 w * za1p * u * za1p * w * zcp * u * zcp +
                                 w * za1p * u * za1p * w * zcp * u * za1p * w, name = "coCC"),
                     # EA
                     mxAlgebra(0, name = "wiEA"),
                     mxAlgebra(0, name = "siEA"),
                     mxAlgebra(zep * u * za1p +
                                 zep * u * zcp * w, name = "paEA"),
                     mxAlgebra(0, name = "spEA"),
                     mxAlgebra(zep * u * za1p * aln +
                                 zep * u * zcp * w, name = "psEA"),
                     mxAlgebra(zep * u * za1p * aln * za1p * u * za1p +
                                 zep * u * za1p * aln * za1p * u * zcp * w +
                                 
                                 zep * u * za1p * w * zcp * u * za1p +
                                 zep * u * za1p * w * zcp * u * zcp * w +
                                 
                                 zep * u * zcp * 1 * zcp * u * za1p +
                                 zep * u * zcp * 1 * zcp * u * zcp * w +
                                 
                                 zep * u * zcp * w * za1p * u * za1p +
                                 zep * u * zcp * w * za1p * u * zcp * w, name = "coEA"),
                     # CA
                     mxAlgebra(w, name = "wiCA"),
                     mxAlgebra(w, name = "siCA"),
                     mxAlgebra(zcp * u * za1p +
                                 zcp * u * zcp * w +
                                 w * za1p * u * za1p +
                                 w * za1p * u * zcp * w, name = "paCA"),
                     mxAlgebra(1 * zcp * u * za1p +
                                 1 * zcp * u * zcp * w +
                                 w * za1p * u * za1p +
                                 w * za1p * u * zcp * w, name = "spCA"),
                     mxAlgebra(zcp * u * za1p * aln +
                                 zcp * u * zcp * w +
                                 w * za1p * u * za1p * aln +
                                 w * za1p * u * zcp * w, name = "psCA"),
                     mxAlgebra(zcp * u * zcp * 1 * zcp * u * za1p +
                                 zcp * u * zcp * 1 * zcp * u * zcp * w +
                                 zcp * u * zcp * w * za1p * u * za1p +
                                 zcp * u * zcp * w * za1p * u * zcp * w +
                                 
                                 zcp * u * za1p * aln * za1p * u * za1p +
                                 zcp * u * za1p * aln * za1p * u * zcp * w +
                                 zcp * u * za1p * w * zcp * u * za1p +
                                 zcp * u * za1p * w * zcp * u * zcp * w +
                                 
                                 w * za1p * u * zcp * 1 * zcp * u * za1p +
                                 w * za1p * u * zcp * 1 * zcp * u * zcp * w +
                                 w * za1p * u * zcp * w * za1p * u * za1p +
                                 w * za1p * u * zcp * w * za1p * u * zcp * w +
                                 
                                 w * za1p * u * za1p * aln * za1p * u * za1p +
                                 w * za1p * u * za1p * aln * za1p * u * zcp * w +
                                 w * za1p * u * za1p * w * zcp * u * za1p +
                                 w * za1p * u * za1p * w * zcp * u * zcp * w, name = "coCA"),
                     # CE
                     mxAlgebra(0, name = "wiCE"),
                     mxAlgebra(0, name = "siCE"),
                     mxAlgebra(zcp * u * zep +
                                 w * za1p * u * zep, name = "paCE"),
                     mxAlgebra(1 * zcp * u * zep +
                                 w * za1p * u * zep, name = "spCE"),
                     mxAlgebra(0, name = "psCE"),
                     mxAlgebra(zcp * u * za1p * aln * za1p * u * zep +
                                 zcp * u * za1p * w * zcp * u * zep +
                                 
                                 zcp * u * zcp * 1 * zcp * u * zep +
                                 zcp * u * zcp * w * za1p * u * zep +
                                 
                                 w * za1p * u * za1p * aln * za1p * u * zep +
                                 w * za1p * u * za1p * w * zcp * u * zep +
                                 
                                 w * za1p * u * zcp * 1 * zcp * u * zep +
                                 w * za1p * u * zcp * w * za1p * u * zep, name = "coCE"),
                     
                     # Now all those correlations need to be positioned correctly
                     # AA matrices
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiAA, siAA),
                                     cbind(siAA, wiAA)), name = "SAAmmMM"),
                     mxAlgebra(rbind(cbind(wiAA, coAA),
                                     cbind(coAA, wiAA)), name = "SAAffMM"),
                     mxAlgebra(rbind(cbind(paAA, psAA),
                                     cbind(psAA, paAA)), name = "SAAfmMM"),
                     mxAlgebra(rbind(cbind(paAA, spAA),
                                     cbind(spAA, paAA)), name = "SAAmfMM"),
                     # Father1 & father2
                     mxAlgebra(rbind(cbind(wiAA, coAA),
                                     cbind(coAA, wiAA)), name = "SAAmmFF"),
                     mxAlgebra(rbind(cbind(wiAA, siAA),
                                     cbind(siAA, wiAA)), name = "SAAffFF"),
                     mxAlgebra(rbind(cbind(paAA, spAA),
                                     cbind(spAA, paAA)), name = "SAAfmFF"),
                     mxAlgebra(rbind(cbind(paAA, psAA),
                                     cbind(psAA, paAA)), name = "SAAmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiAA, spAA),
                                     cbind(psAA, wiAA)), name = "SAAmmMF"),
                     mxAlgebra(rbind(cbind(wiAA, psAA),
                                     cbind(spAA, wiAA)), name = "SAAffMF"),
                     mxAlgebra(rbind(cbind(paAA, coAA),
                                     cbind(siAA, paAA)), name = "SAAfmMF"),
                     mxAlgebra(rbind(cbind(paAA, siAA),
                                     cbind(coAA, paAA)), name = "SAAmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiAA, psAA),
                                     cbind(spAA, wiAA)), name = "SAAmmFM"),
                     mxAlgebra(rbind(cbind(wiAA, spAA),
                                     cbind(psAA, wiAA)), name = "SAAffFM"),
                     mxAlgebra(rbind(cbind(paAA, siAA),
                                     cbind(coAA, paAA)), name = "SAAfmFM"),
                     mxAlgebra(rbind(cbind(paAA, coAA),
                                     cbind(siAA, paAA)), name = "SAAmfFM"),
                     
                     # EE
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiEE, siEE),
                                     cbind(siEE, wiEE)), name = "SEEmmMM"),
                     mxAlgebra(rbind(cbind(wiEE, coEE),
                                     cbind(coEE, wiEE)), name = "SEEffMM"),
                     mxAlgebra(rbind(cbind(paEE, psEE),
                                     cbind(psEE, paEE)), name = "SEEfmMM"),
                     mxAlgebra(rbind(cbind(paEE, spEE),
                                     cbind(spEE, paEE)), name = "SEEmfMM"),
                     # Father1 & father2
                     mxAlgebra(rbind(cbind(wiEE, coEE),
                                     cbind(coEE, wiEE)), name = "SEEmmFF"),
                     mxAlgebra(rbind(cbind(wiEE, siEE),
                                     cbind(siEE, wiEE)), name = "SEEffFF"),
                     mxAlgebra(rbind(cbind(paEE, spEE),
                                     cbind(spEE, paEE)), name = "SEEfmFF"),
                     mxAlgebra(rbind(cbind(paEE, psEE),
                                     cbind(psEE, paEE)), name = "SEEmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiEE, spEE),
                                     cbind(psEE, wiEE)), name = "SEEmmMF"),
                     mxAlgebra(rbind(cbind(wiEE, psEE),
                                     cbind(spEE, wiEE)), name = "SEEffMF"),
                     mxAlgebra(rbind(cbind(paEE, coEE),
                                     cbind(siEE, paEE)), name = "SEEfmMF"),
                     mxAlgebra(rbind(cbind(paEE, siEE),
                                     cbind(coEE, paEE)), name = "SEEmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiEE, psEE),
                                     cbind(spEE, wiEE)), name = "SEEmmFM"),
                     mxAlgebra(rbind(cbind(wiEE, spEE),
                                     cbind(psEE, wiEE)), name = "SEEffFM"),
                     mxAlgebra(rbind(cbind(paEE, siEE),
                                     cbind(coEE, paEE)), name = "SEEfmFM"),
                     mxAlgebra(rbind(cbind(paEE, coEE),
                                     cbind(siEE, paEE)), name = "SEEmfFM"),
                     
                     # CC
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiCC, siCC),
                                     cbind(siCC, wiCC)), name = "SCCmmMM"),
                     mxAlgebra(rbind(cbind(wiCC, coCC),
                                     cbind(coCC, wiCC)), name = "SCCffMM"),
                     mxAlgebra(rbind(cbind(paCC, psCC),
                                     cbind(psCC, paCC)), name = "SCCfmMM"),
                     mxAlgebra(rbind(cbind(paCC, spCC),
                                     cbind(spCC, paCC)), name = "SCCmfMM"),
                     # Father1 & father2
                     mxAlgebra(rbind(cbind(wiCC, coCC),
                                     cbind(coCC, wiCC)), name = "SCCmmFF"),
                     mxAlgebra(rbind(cbind(wiCC, siCC),
                                     cbind(siCC, wiCC)), name = "SCCffFF"),
                     mxAlgebra(rbind(cbind(paCC, spCC),
                                     cbind(spCC, paCC)), name = "SCCfmFF"),
                     mxAlgebra(rbind(cbind(paCC, psCC),
                                     cbind(psCC, paCC)), name = "SCCmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiCC, spCC),
                                     cbind(psCC, wiCC)), name = "SCCmmMF"),
                     mxAlgebra(rbind(cbind(wiCC, psCC),
                                     cbind(spCC, wiCC)), name = "SCCffMF"),
                     mxAlgebra(rbind(cbind(paCC, coCC),
                                     cbind(siCC, paCC)), name = "SCCfmMF"),
                     mxAlgebra(rbind(cbind(paCC, siCC),
                                     cbind(coCC, paCC)), name = "SCCmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiCC, psCC),
                                     cbind(spCC, wiCC)), name = "SCCmmFM"),
                     mxAlgebra(rbind(cbind(wiCC, spCC),
                                     cbind(psCC, wiCC)), name = "SCCffFM"),
                     mxAlgebra(rbind(cbind(paCC, siCC),
                                     cbind(coCC, paCC)), name = "SCCfmFM"),
                     mxAlgebra(rbind(cbind(paCC, coCC),
                                     cbind(siCC, paCC)), name = "SCCmfFM"),
                     
                     # EA
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiEA, siEA),
                                     cbind(siEA, wiEA)), name = "SEAmmMM"),
                     mxAlgebra(rbind(cbind(wiEA, coEA),
                                     cbind(coEA, wiEA)), name = "SEAffMM"),
                     mxAlgebra(rbind(cbind(paEA, psEA),
                                     cbind(psEA, paEA)), name = "SEAfmMM"),
                     mxAlgebra(rbind(cbind(paEA, spEA),
                                     cbind(spEA, paEA)), name = "SEAmfMM"),
                     # Father1 & father2
                     mxAlgebra(rbind(cbind(wiEA, coEA),
                                     cbind(coEA, wiEA)), name = "SEAmmFF"),
                     mxAlgebra(rbind(cbind(wiEA, siEA),
                                     cbind(siEA, wiEA)), name = "SEAffFF"),
                     mxAlgebra(rbind(cbind(paEA, spEA),
                                     cbind(spEA, paEA)), name = "SEAfmFF"),
                     mxAlgebra(rbind(cbind(paEA, psEA),
                                     cbind(psEA, paEA)), name = "SEAmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiEA, spEA),
                                     cbind(psEA, wiEA)), name = "SEAmmMF"),
                     mxAlgebra(rbind(cbind(wiEA, psEA),
                                     cbind(spEA, wiEA)), name = "SEAffMF"),
                     mxAlgebra(rbind(cbind(paEA, coEA),
                                     cbind(siEA, paEA)), name = "SEAfmMF"),
                     mxAlgebra(rbind(cbind(paEA, siEA),
                                     cbind(coEA, paEA)), name = "SEAmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiEA, psEA),
                                     cbind(spEA, wiEA)), name = "SEAmmFM"),
                     mxAlgebra(rbind(cbind(wiEA, spEA),
                                     cbind(psEA, wiEA)), name = "SEAffFM"),
                     mxAlgebra(rbind(cbind(paEA, siEA),
                                     cbind(coEA, paEA)), name = "SEAfmFM"),
                     mxAlgebra(rbind(cbind(paEA, coEA),
                                     cbind(siEA, paEA)), name = "SEAmfFM"),
                     
                     # CA
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiCA, siCA),
                                     cbind(siCA, wiCA)), name = "SCAmmMM"),
                     mxAlgebra(rbind(cbind(wiCA, coCA),
                                     cbind(coCA, wiCA)), name = "SCAffMM"),
                     mxAlgebra(rbind(cbind(paCA, psCA),
                                     cbind(psCA, paCA)), name = "SCAfmMM"),
                     mxAlgebra(rbind(cbind(paCA, spCA),
                                     cbind(spCA, paCA)), name = "SCAmfMM"),
                     # father1 & father2
                     mxAlgebra(rbind(cbind(wiCA, coCA),
                                     cbind(coCA, wiCA)), name = "SCAmmFF"),
                     mxAlgebra(rbind(cbind(wiCA, siCA),
                                     cbind(siCA, wiCA)), name = "SCAffFF"),
                     mxAlgebra(rbind(cbind(paCA, spCA),
                                     cbind(spCA, paCA)), name = "SCAfmFF"),
                     mxAlgebra(rbind(cbind(paCA, psCA),
                                     cbind(psCA, paCA)), name = "SCAmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiCA, spCA),
                                     cbind(psCA, wiCA)), name = "SCAmmMF"),
                     mxAlgebra(rbind(cbind(wiCA, psCA),
                                     cbind(spCA, wiCA)), name = "SCAffMF"),
                     mxAlgebra(rbind(cbind(paCA, coCA),
                                     cbind(siCA, paCA)), name = "SCAfmMF"),
                     mxAlgebra(rbind(cbind(paCA, siCA),
                                     cbind(coCA, paCA)), name = "SCAmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiCA, psCA),
                                     cbind(spCA, wiCA)), name = "SCAmmFM"),
                     mxAlgebra(rbind(cbind(wiCA, spCA),
                                     cbind(psCA, wiCA)), name = "SCAffFM"),
                     mxAlgebra(rbind(cbind(paCA, siCA),
                                     cbind(coCA, paCA)), name = "SCAfmFM"),
                     mxAlgebra(rbind(cbind(paCA, coCA),
                                     cbind(siCA, paCA)), name = "SCAmfFM"),
                     
                     # CE
                     # ------------------------------------------------------------
                     # Mother1 & mother2
                     mxAlgebra(rbind(cbind(wiCE, siCE),
                                     cbind(siCE, wiCE)), name = "SCEmmMM"),
                     mxAlgebra(rbind(cbind(wiCE, coCE),
                                     cbind(coCE, wiCE)), name = "SCEffMM"),
                     mxAlgebra(rbind(cbind(paCE, psCE),
                                     cbind(psCE, paCE)), name = "SCEfmMM"),
                     mxAlgebra(rbind(cbind(paCE, spCE),
                                     cbind(spCE, paCE)), name = "SCEmfMM"),
                     # father1 & father2
                     mxAlgebra(rbind(cbind(wiCE, coCE),
                                     cbind(coCE, wiCE)), name = "SCEmmFF"),
                     mxAlgebra(rbind(cbind(wiCE, siCE),
                                     cbind(siCE, wiCE)), name = "SCEffFF"),
                     mxAlgebra(rbind(cbind(paCE, spCE),
                                     cbind(spCE, paCE)), name = "SCEfmFF"),
                     mxAlgebra(rbind(cbind(paCE, psCE),
                                     cbind(psCE, paCE)), name = "SCEmfFF"),
                     # Mother1 & father2
                     mxAlgebra(rbind(cbind(wiCE, spCE),
                                     cbind(psCE, wiCE)), name = "SCEmmMF"),
                     mxAlgebra(rbind(cbind(wiCE, psCE),
                                     cbind(spCE, wiCE)), name = "SCEffMF"),
                     mxAlgebra(rbind(cbind(paCE, coCE),
                                     cbind(siCE, paCE)), name = "SCEfmMF"),
                     mxAlgebra(rbind(cbind(paCE, siCE),
                                     cbind(coCE, paCE)), name = "SCEmfMF"),
                     # Father1 & mother2
                     mxAlgebra(rbind(cbind(wiCE, psCE),
                                     cbind(spCE, wiCE)), name = "SCEmmFM"),
                     mxAlgebra(rbind(cbind(wiCE, spCE),
                                     cbind(psCE, wiCE)), name = "SCEffFM"),
                     mxAlgebra(rbind(cbind(paCE, siCE),
                                     cbind(coCE, paCE)), name = "SCEfmFM"),
                     mxAlgebra(rbind(cbind(paCE, coCE),
                                     cbind(siCE, paCE)), name = "SCEmfFM"),
                     
                     # Build S matrix
                     # ------------------------------------------------------------
                     mxAlgebra(SAAmmMM * data.mm + SAAmmFF * data.ff + SAAmmMF * data.mf + SAAmmFM * data.fm, name = "SAAmm"),
                     mxAlgebra(SAAfmMM * data.mm + SAAfmFF * data.ff + SAAfmMF * data.mf + SAAfmFM * data.fm, name = "SAAfm"),
                     mxAlgebra(SEAmmMM * data.mm + SEAmmFF * data.ff + SEAmmMF * data.mf + SEAmmFM * data.fm, name = "SEAmm"),
                     mxAlgebra(SEAfmMM * data.mm + SEAfmFF * data.ff + SEAfmMF * data.mf + SEAfmFM * data.fm, name = "SEAfm"),
                     mxAlgebra(SCAmmMM * data.mm + SCAmmFF * data.ff + SCAmmMF * data.mf + SCAmmFM * data.fm, name = "SCAmm"),
                     mxAlgebra(SCAfmMM * data.mm + SCAfmFF * data.ff + SCAfmMF * data.mf + SCAfmFM * data.fm, name = "SCAfm"),
                     
                     mxAlgebra(SAAffMM * data.mm + SAAffFF * data.ff + SAAffMF * data.mf + SAAffFM * data.fm, name = "SAAff"),
                     mxAlgebra(SEAmfMM * data.mm + SEAmfFF * data.ff + SEAmfMF * data.mf + SEAmfFM * data.fm, name = "SEAmf"),
                     mxAlgebra(SEAffMM * data.mm + SEAffFF * data.ff + SEAffMF * data.mf + SEAffFM * data.fm, name = "SEAff"),
                     mxAlgebra(SCAmfMM * data.mm + SCAmfFF * data.ff + SCAmfMF * data.mf + SCAmfFM * data.fm, name = "SCAmf"),
                     mxAlgebra(SCAffMM * data.mm + SCAffFF * data.ff + SCAffMF * data.mf + SCAffFM * data.fm, name = "SCAff"),
                     
                     mxAlgebra(SEEmmMM * data.mm + SEEmmFF * data.ff + SEEmmMF * data.mf + SEEmmFM * data.fm, name = "SEEmm"),
                     mxAlgebra(SEEfmMM * data.mm + SEEfmFF * data.ff + SEEfmMF * data.mf + SEEfmFM * data.fm, name = "SEEfm"),
                     mxAlgebra(SCEmmMM * data.mm + SCEmmFF * data.ff + SCEmmMF * data.mf + SCEmmFM * data.fm, name = "SCEmm"),
                     mxAlgebra(SCEfmMM * data.mm + SCEfmFF * data.ff + SCEfmMF * data.mf + SCEfmFM * data.fm, name = "SCEfm"),
                     
                     mxAlgebra(SEEffMM * data.mm + SEEffFF * data.ff + SEEffMF * data.mf + SEEffFM * data.fm, name = "SEEff"),
                     mxAlgebra(SCEmfMM * data.mm + SCEmfFF * data.ff + SCEmfMF * data.mf + SCEmfFM * data.fm, name = "SCEmf"),
                     mxAlgebra(SCEffMM * data.mm + SCEffFF * data.ff + SCEffMF * data.mf + SCEffFM * data.fm, name = "SCEff"),
                     
                     mxAlgebra(SCCmmMM * data.mm + SCCmmFF * data.ff + SCCmmMF * data.mf + SCCmmFM * data.fm, name = "SCCmm"),
                     mxAlgebra(SCCfmMM * data.mm + SCCfmFF * data.ff + SCCfmMF * data.mf + SCCfmFM * data.fm, name = "SCCfm"),
                     
                     mxAlgebra(SCCffMM * data.mm + SCCffFF * data.ff + SCCffMF * data.mf + SCCffFM * data.fm, name = "SCCff"),
                     
                     # Y_mm, Y_ff, A1_mm, A1_ff, E_mm, E_ff, C_mm, C_ff
                     mxAlgebra(rbind(cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, SAAmm, t(SAAfm), t(SEAmm), t(SEAfm), t(SCAmm), t(SCAfm)),
                                     cbind(Z, Z, SAAfm, SAAff, t(SEAmf), t(SEAff), t(SCAmf), t(SCAff)),
                                     cbind(Z, Z, SEAmm, SEAmf, SEEmm, t(SEEfm), t(SCEmm), t(SCEfm)),
                                     cbind(Z, Z, SEAfm, SEAff, SEEfm, SEEff, t(SCEmf), t(SCEff)),
                                     cbind(Z, Z, SCAmm, SCAmf, SCEmm, SCEmf, SCCmm, t(SCCfm)),
                                     cbind(Z, Z, SCAfm, SCAff, SCEfm, SCEff, SCCfm, SCCff)), name = "Spp"),
                     
                     # A matrix
                     mxAlgebra(rbind(cbind(Z, Z, Aa1p, Z, Aep, Z, Acp, Z),
                                     cbind(Z, Z, Z, Aa1p, Z, Aep, Z, Acp),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z)), name = "App"),
                     # V
                     mxAlgebra(solve(Ipp - App) %*% Spp %*% t(solve(Ipp - App)), name = "Vpp"),
                     
                     # Means
                     mxAlgebra(rbind(Um, Uf, Z1, Z1, Z1, Z1, Z1, Z1), name = "Upp"),
                     mxAlgebra(solve(Ipp - App) %*% Upp, name = "Mpp"),
                     
                     # Offspring model
                     # -----------------------------------------------------------------------------------
                     # A - op
                     mxAlgebra(rbind(cbind(Abm, Abf, Z, Z, Z, Z, Z, Z),
                                     cbind(Abm, Abf, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Atr, Atr, Z, Z, Z, Z),
                                     cbind(Z, Z, Atr, Atr, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z)), name = "Aop"),
                     # Sr - oo
                     mxAlgebra(rbind(cbind(1 - 0.5 * (1 + paAA), 0),
                                     cbind(0, 1 - 0.5 * (1 + paAA))), name = "SA1oo"),
                     mxAlgebra(rbind(cbind(1, 0.5^2 * al),
                                     cbind(0.5^2 * al, 1)), name = "SA2ooB"),
                     mxAlgebra(rbind(cbind(0.5, 0.5^2 * al),
                                     cbind(0.5^2 * al, 0.5)), name = "SA2ooW"),
                     mxAlgebra(rbind(cbind(1, 0),
                                     cbind(0, 1)), name = "SEoo"),
                     # C
                     mxAlgebra(exp(-0.5 * (abs(data.cohort_o11 - data.cohort_o21) / rho)^2), name = "rCo1"),
                     mxAlgebra(exp(-0.5 * (abs(data.cohort_o12 - data.cohort_o22) / rho)^2), name = "rCo2"),
                     mxAlgebra(rbind(cbind(1, 0),
                                      cbind(0, 1)), name = "SCooB"),
                     mxAlgebra(rbind(cbind(rCo1, 0),
                                     cbind(0, rCo2)), name = "SCooW"),
                     #mxAlgebra(rbind(cbind(1, 0),
                    #                 cbind(0, 1)), name = "SCooW"),
                     # Y_oo, A1_oo, A2_oo, E_oo
                     mxAlgebra(rbind(cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, SA1oo, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, SA1oo, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, SA2ooB, SA2ooW, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, SA2ooW, SA2ooB, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, SEoo, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, SEoo, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, SCooB, SCooW),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, SCooW, SCooB)), name = "Sroo"),
                     # S - oo
                     mxAlgebra(Aop %*% Vpp %*% t(Aop) + Sroo, name = "Soo"),
                     # A - oo
                     mxAlgebra(rbind(cbind(Z, Z, Aa1o, Z, Aa2o, Z, Aeo, Z, Aco, Z),
                                     cbind(Z, Z, Z, Aa1o, Z, Aa2o, Z, Aeo, Z, Aco),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, Z, Z, Z, Z, Z, Z, Z, Z, Z)), name = "Aoo"),
                     # V - oo
                     mxAlgebra(solve(Ioo - Aoo) %*% Soo %*% t(solve(Ioo - Aoo)), name = "Voo"),
                     
                     # Means
                     mxAlgebra(rbind(Uo + Bfemaleo %*% sex_o1, Uo + Bfemaleo %*% sex_o2, Z1, Z1, Z1, Z1, Z1, Z1, Z1, Z1), name = "Uoo"),
                     mxAlgebra(solve(Ioo - Aoo) %*% Uoo + Aop %*% Mpp, name = "Moo"),
                     
                     # Manifests
                     # ---------------------------------------------------------------------------------
                     # Filters
                     mxAlgebra(rbind(cbind(I, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, I, Z, Z, Z, Z, Z, Z)), name = "Fpp"),
                     mxAlgebra(rbind(cbind(I, Z, Z, Z, Z, Z, Z, Z, Z, Z),
                                     cbind(Z, I, Z, Z, Z, Z, Z, Z, Z, Z)), name = "Foo"),
                     # Manifest covariance
                     mxAlgebra(Foo %*% Voo %*% t(Foo), name = "Coo"),
                     mxAlgebra(Fpp %*% Vpp %*% t(Fpp), name = "Cpp"),
                     mxAlgebra(Fpp %*% Vpp %*% t(Aop) %*% t(solve(Ioo - Aoo)) %*% t(Foo), name = "Cpo"),
                     mxAlgebra(Foo %*% solve(Ioo - Aoo) %*% Aop %*% t(Vpp) %*% t(Fpp), name = "Cop"),
                     mxAlgebra(rbind(cbind(Cpp, Cpo),
                                     cbind(Cop, Coo)), name = "C"),
                     # Manifest means
                     mxAlgebra(Fpp %*% Mpp, name = "Mmpp"),
                     mxAlgebra(Foo %*% Moo, name = "Mmoo"),
                     mxAlgebra(t(rbind(Mmpp, Mmoo)), name = "M"),
                     # Confidence intervals - THL
                     mxCI(c("a1p","a1o","a2o","cp","ep","eo","co","p","m","d","w","um","uf","uo","bfemaleo","rho"), interval=0.95, type="both"),
                     # Expectation
                     mxExpectationNormal("C", "M", var_nms[1:8]),
                     # Objective
                     mxFitFunctionML())
  mod_full_ci <- mxModel(mod_full, mxCI(c("m","p","v2p","d"), type="both",interval=0.95))
  mxOption(mod_full_ci, "Calculate Hessian","Yes")
  fit_full1 = mxRun(mod_full_ci)
  result <- fit_full1
  # if (tryhard>0) {
  #   fit_full12 = mxTryHardWideSearch(mod_full_ci, extraTries = tryhard)
  #   result <- fit_full12
  # }
  return(result)
}

#summary(fit_full12)
#round(coef(fit_full12), 3)
# # Only genetic transmission
# mod_gen = omxSetParameters(mxModel(mod_full, "Equate_rae_gen", remove = T), c("m", "p", "w"), F, 0, name = "mod_gen")
# fit_gen1 = mxTryHardWideSearch(mod_gen, extraTries = 6)
# fit_gen12 = mxTryHard(fit_gen1)
# summary(fit_gen12)
# round(coef(fit_gen12), 3)
# 
# # Only environmental transmission
# mod_env = omxSetParameters(mxModel(mod_full, "Equate_rae_gen", remove = T), c("a1o", "w"), F, 0, name = "mod_env")
# fit_env1 = mxTryHardWideSearch(mod_env, extraTries = 5)
# fit_env12 = mxTryHard(fit_env1)
# summary(fit_env12)
# round(coef(fit_env12), 3)
# 
# # Compare models
# mxCompare(fit_full12, list(fit_gen12, fit_env12))

# res =list(fit_full = fit_full12,
#           fit_gen = fit_gen12,
#           fit_env = fit_env12) 

#save(res, file = "N:/data/durable/projects/intertranseduc/results/fit_ACEpACEc_rCA_STD.RData")


