function [VaR, ETL] = cornishFisher(returns, conf, inv0)
%Function to Calculate VaR and ETL By Cornish-Fisher Approximation.

%Takes a Vector of PnL and a confidence level (i.e 0.95 or 0.99)
%and calculates VaR & ETLusing a 4th order CF approximation.

%INPUTS:
%       returns: (nx1): Vector of Portfolio profit & loss
%       conf: (scalar); Confidence Level for VaR & ETL (i.e 0.95, 0.99)
%       inv0: (scalar); Investment or Portfolio Size

%OUTPUTS:
%       VaR (scalar): Value at Risk under CF approximation
%       ETL (scalar): Expected Tail Loss under CF approximation

%***************************************************************************************

%Calculate First Four Moments of Empirical Distribution:
mu = mean(returns);
sigma = std(returns);
skew = skewness(returns);
exKurtosis = kurtosis(returns) - 3; 


%ath quantile of the Standard Normal:
a = 1-conf;
z = norminv(a, 0,1);

%Fourth order CF expansion:
tau = skew/6;
xa = z + tau*(z^2 - 1) + (exKurtosis/24)*z*(z^2 - 3) - (tau^2)*z*(2*z^2 - 5);

%Cornish-Fisher VaR:
VaR = -(xa*sigma + mu)*inv0;
%ETL: (as in Alexander pg 195):
x2 = (1/a)*normpdf(z, 0,1);
f = x2 + tau*(x2^2 - 1) + (exKurtosis/24)*x2*(x2^2 - 3) - (tau^2)*x2*(2*x2^2 - 5);

ETL = (f*sigma - mu)*inv0;

end 


