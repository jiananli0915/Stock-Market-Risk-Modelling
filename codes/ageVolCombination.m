function [VaR, ETL] = ageVolCombination(returns, windowSize, conf, inv0, forecastVec, lambda, trainingSize, model, P, Q, nu)
%ageVolCombination: Function to calculate a rolling Value-At-Risk By
%Weighting Losses By Volatility And Time of Occurence.

%This function calculates a rolling VaR by weighting a sample of returns by
%volatility and weighting these returns by age.

%The method used to forecast Volatility is GARCH. 
%The method used to weight time is EWMA.

%INPUTS:
%       returns: (nx1); Vector of returns.
%       windowSize: (scalar); Size of estimation window for calculating
%                               VaRs
%       forecastVec: (kx1); Vector of Rolling Volatility Forecasts (more
%                       efficient than calculating in-function
%       lambda: Decay factor for EWMA
%       trainingSize: (scalar); Number of samples used to calibrate the
%                               GARCH model
%       model: 
%                   MUST BE IN [0,1,2]--> encoded data to specify which GARCH model to use.
%                   0 = GARCH(P,Q)
%                   1 = EGARCH(P,Q)
%                   2 = GJR(P,Q)
%                   3 = GARCH(P,Q) (Under T with nu Degrees of Freedom
%       lambda: (scalar); Weighting Factor
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Confidence Level for VaR (i.e 0.95, 0.99)

%OUTPUTS:
%       VaR (nEstimates); Rolling Value-At-Risk
%       estimates at confidence level 'conf'
%       ETL (nEstimates); Rolling estimates at confidence level 'conf'
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
   
    [VaR(i), ETL(i)] = ageWeightedEWMA(weightedReturns, lambda, inv0, conf);

end 


