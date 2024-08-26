# EPD under Progressive Type-II Censuring with Binomial Removals

This is a sample code from the main research project on Exponential Power Distribution. It is a sample demonstration on how to estimate the parameters (alpha and lambda) of an EPD under Progressive Type-II Censuring with Binomial Removals using Metropolis Hastings Algorithm.

We will be performing the following operations in our code:-

    - Generate a progressive Type-II censored sample with binomial removals from the EPD.
    - Set the values of the hyperparameters in such a manner that the prior distributions non-informative. 
    - Then we will generate sample observations from posterior distributions using Markov Chain Monte Carlo methodology for alpha and lambda - (Metropolis Hastings algorithm in this case).
    - Generarte sample observations for 'p' using standard generation method for first kind beta distribution (i.e. we will use ziggurat algorithm - rbeta() function  in R)
    - We will obtain Bayes Estimates under SELF-loss function.
    - For Bayesian Credible Intervals (BCI), we are computing Equal-tailed Credible Interval (ETCI) and Highest Posterior Density Credible Interval (HPDCI) and we are also calculating their respective widths and Coverage Probabilities (CP)
    - Also, we will be calcuating Average Absolute Bias (AB) and Mean Square Error (MSE) by repeating the steps above.

Note:- 
(1) We will be using mcmc library to perform M-H algorithm and we will be including burn-in as well and maintaining an acceptance sampling rate of 20%
(2) With gamma priors, posterior distribution of alpha and lambda are still not in standard form (cannot be solved analytically due to multi-dimensional integration), hence we will use Metropolis-Hastings algorithm to generate sammple observations from posterior distributions. 
(3) With beta priors, posterior distribution of 'p'(conditional distribution of 'p' given alpha,lambda and data) follows beta first kind distribution (beta is conjugate prior for binomial), hence we can use standard random generation method. Also, here we will choose Jeffreys Prior as our non-informative prior.
