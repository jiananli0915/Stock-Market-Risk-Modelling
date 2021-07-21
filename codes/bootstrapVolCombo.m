function [VaR, ETL] = bootstrapVolCombo(returns, windowSize, conf, inv0, forecastVec, nboots, trainingSize, model, P, Q, nu)
%bootstrapVolCombo: Function to calculate a rolling Value-At-Risk By
%Weighting Losses By Volatility Combined with Bootstrap Aggregation


%This function calculates a rolling VaR by weighting a sample of returns by
%volatility resampling from each estimation window

%The method used to forecast Volatility is GARCH. 

%General Algorithm:
%       -define estimation window etc
%       -index appropriately to get the correct sample returns
%               -can be checked at the end, i.e have we estimated up to 1st
%               Jan 2020? If all estimates FINISH here we can just look
%               backward and take the appropriate sample. 
%       -estimate a garch model for this period & get one-step ahead
%       forecast
%       -using this forecast, calculate volatility weighted returns 
%       -randomly sample (with replacement) from these returns
%       -shift forward one period, repeat 
%       -end 

%INPUTS:
%       returns: (nx1); Vector of returns.
%       windowSize: (scalar); Size of estimation window for calculating
%                               VaRs
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Confidence Level for VaR (i.e 0.95, 0.99)
%       forecastVec: (kx1); Vector of Rolling Volatility Forecasts (more
%                       efficient than calculating in-function
%       nboots: Number of Bootstrap samples
%       trainingSize: (scalar); Number of samples used to calibrate the
%                               GARCH model
%       model: 
%                   MUST BE IN [0,1,2]--> encoded data to specify which GARCH model to use.
%                   0 = GARCH(P,Q)
%                   1 = EGARCH(P,Q)
%                   2 = GJR(P,Q)
%                   3 = GARCH(P,Q) (Under T with nu Degrees of Freedom
%       lambda: (scalar); Weighting Factor


%OUTPUTS:
%       VaR (nEstimates); Rolling VaR Estimate
%***************************************************************************************

%More efficient to pass a pre-calculated forecastVec, but if we want to use
%a different model, pass extra parameters:
if nargin > 6
    forecastVec = garchForecast(returns, model, nu, P, Q,  trainingSize);
end 

%Returns left over after forecasting:
returnVec = returns(end-length(forecastVec)+2:end);

%Number of VaR Estimates:
nEstimates = length(returnVec) - windowSize + 1;

%Initialize Vectors:
VaR = zeros(nEstimates, 1);
ETL = zeros(nEstimates, 1);

for i = 1:nEstimates 
    
   %pull a sample of returns that is of length 'windowSize':
   estimationReturns = returnVec(i: i+ windowSize - 1);
   estimationForecast = forecastVec(i: i+windowSize -1);
   
   %scale by the 'most recent' volatility forecast:
   sigmaT= forecastVec(i + windowSize);
   
   %Weight estimation window returns by volatility:
   weightedReturns = estimationReturns.*(sigmaT./estimationForecast);
   
    [VaR(i)] = bootstrapHistoricalVaR(weightedReturns, inv0, conf, nboots, 0);
    [ETL(i)] = bootstrapHistoricalVaR(weightedReturns, inv0, conf, nboots, 1);
  

end 

