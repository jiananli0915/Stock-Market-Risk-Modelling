function [VaR, ETL, weightedReturns] = volWeightedVaR(returns, forecastVec, conf, inv0)

%Returns left over after forecasting:
returnVec = returns(end-length(forecastVec)+2:end);

%Number of returns we are weighting:
n = length(returnVec);

%Initialize Vectors:
VaR = zeros(n, 1);
ETL = zeros(n,1);

%Most recent volatility estimate
sigmaT = forecastVec(end);

%weighting returns 
weightedReturns = returnVec.*(sigmaT./forecastVec(1:end-1));

%VaR & ETL:
VaR= historicalVaR(weightedReturns, conf, inv0);
ETL = historicalETL(weightedReturns, conf, inv0 );

end

