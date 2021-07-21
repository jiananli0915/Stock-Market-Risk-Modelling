function [indicator, realizedVolatility] = exceedanceIndicator(returns, VaR)
%exceedanceIndicator: Function to Calculate the Number of Times a Rolling
%Value-At-Risk Estimate is Breached.

%INPUTS:
%           returns: (nx1); PnL Vector
%           VaR: (kx1); Rolling Value-At-Risk Estimates

%           These should both be in the same units***

%OUTPUTS:
%           indictator: (kx1); Logical Vector Taking a Value of 1 if An
%           Exceedance Occurs, 0 otherwise

%*************************************************************************


sampleReturns = returns(end - length(VaR) + 1:end)
realizedVolatility = sqrt(returns.^2);

indicator = (sampleReturns > VaR)
end

