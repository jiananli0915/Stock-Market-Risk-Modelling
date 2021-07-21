function [VaR] = rollingExtremeVaR(returns, conf, inv0, windowSize)


nEstimates = length(returns) - windowSize;
VaR = zeros(nEstimates, 1);

for j = 1:nEstimates 
    t = windowSize + j;
    estWindow = returns(t-499:t);
 
    VaR(j) =evVaR(estWindow,inv0, conf);
end

