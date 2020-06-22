# dynquant
Pacakge dynquant implements dynamic quantile models allowing for covariates to be sampled at mixed frequency. Implementation contains the following models:

1. CAViaR (different specifications, including MIDAS data), see [1]. 
2. MVMQCAViaR (including MIDAS data), see [3].
3. ARDL-MIDAS quantile and MIDAS quantile regression models, see *midasml* package, and [2, 4, 5].

##
The main functions are ```fit_caviar``` (CAViaR) & ```fit_mvmqcaviar``` (MVMQCAViaR) - both have options to fit MIDAS data. 


# References

[1] Engle, R. F., & Manganelli, S. (2004). CAViaR: Conditional autoregressive value at risk by regression quantiles. Journal of Business & Economic Statistics, 22(4), 367-381. https://doi.org/10.1198/073500104000000370

[2] Ghysels, E. (2014). Conditional skewness with quantile regression models: SoFiE Presidential Address and a tribute to Hal White. Journal of Financial Econometrics, 12(4), 620-644. https://doi.org/10.1093/jjfinec/nbu021

[3] White, H., Kim, T. H., & Manganelli, S. (2015). VAR for VaR: Measuring tail dependence using multivariate regression quantiles. Journal of Econometrics, 187(1), 169-188. https://doi.org/10.1016/j.jeconom.2015.02.004

[4] Ghysels, E., Plazzi, A., & Valkanov, R. (2016). Why invest in emerging markets? The role of conditional return asymmetry. The Journal of Finance, 71(5), 2145-2192. https://doi.org/10.1111/jofi.12420

[5] Ghysels, E., Iania, L., Striaukas, J. (2018). Quantile-based Inflation Risk Models https://www.nbb.be/doc/ts/publications/wp/wp349en.pdf


