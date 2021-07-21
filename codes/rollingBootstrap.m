function [VaR, ETL] = rollingBootstrap(returns, conf, inv0, windowSize, nboots)
 

%Number of VaR Estimates:
nEstimates = length(returns) - windowSize + 1;

%Initialize Output Vectors:
VaR = zeros(nEstimates, 1);
ETL = zeros(nEstimates, 1);


%Loop to Calculate Rolling VaR Estimate:
for i = 1:nEstimates 
    t = windowSize + i - 1; 

    estWindow = returns(i:t);

    VaR(i) = bootstrapHistoricalVaR(estWindow, inv0, conf, nboots, 0);
    ETL(i) = bootstrapHistoricalVaR(estWindow, inv0, conf, nboots, 1);

end 


