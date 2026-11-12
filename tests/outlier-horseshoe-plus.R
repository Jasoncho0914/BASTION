library(BASTION)
ns<-asNamespace('BASTION')
init<-get('init_Outlier',ns);update<-get('fit_Outlier',ns)
scales<-get('t_sampleEvolZeta_ps',ns)
obs<-list(sigma_e=2,sigma_et=rep(2,12))
set.seed(91);p<-init(c(12,-10,9,-8,rep(0,8)),obs)
stopifnot(p$Td==12,length(p$s_mu)==12,all(p$s_mu[1:4]!=0),
  all(c('eta_aux','eta_t2') %in% names(p$s_evolParams)))
# Check every conditional against an independently written seeded sweep.
z<-seq(-1,1,length.out=12);q<-p$s_evolParams
ig<-extraDistr::rinvgamma
set.seed(92)
expected<-q
expected$lambda_t2<-ig(12,1,1/q$v+z^2/(2*q$tau_t2))
expected$v<-ig(12,1,1/expected$lambda_t2+1/q$eta_t2)
expected$eta_t2<-ig(12,1,1/q$eta_aux+1/expected$v)
expected$eta_aux<-ig(12,1,1+1/expected$eta_t2)
expected$tau_t2<-ig(1,6.5,1/q$xi+sum(z^2/expected$lambda_t2)/2)
expected$xi<-ig(1,1,1+1/expected$tau_t2)
expected$sigma_wt<-sqrt(expected$tau_t2*expected$lambda_t2)
set.seed(92);actual<-scales(z,12,q)
stopifnot(isTRUE(all.equal(actual,expected,tolerance=1e-14)))
for(i in 1:100) p<-update(c(12,-10,9,-8,rep(0,8)),p,obs)
stopifnot(p$Td==12,all(is.finite(p$s_mu)),
  isTRUE(all.equal(p$n_squared_sum,sum((p$s_mu/p$s_evolParams$sigma_wt)^2))))
# Invalid state fails explicitly rather than being truncated.
bad<-q;bad$v[1]<-0
stopifnot(inherits(tryCatch(scales(z,12,bad),error=identity),'error'))
# Validate the Gaussian conditional by repeated simulation.
set.seed(94)
draw<-get('sampleOutlier',ns)(rep(3,40000),4,9,40000)
stopifnot(abs(mean(draw)-27/13)<.04,abs(var(draw)-36/13)<.08)
# Constant observation variance counts all T outlier coefficients.
variance_update<-get('fit_sigmaE_0_m',ns)
set.seed(95);reference<-sqrt(ig(1,12,(sum(p$s_mu^2)+p$n_squared_sum)/2))
set.seed(95);actual<-variance_update(rep(0,12),list(p),12)$sigma_e
stopifnot(identical(actual,reference))
cat('Horseshoe+ outlier checks passed.\n')
