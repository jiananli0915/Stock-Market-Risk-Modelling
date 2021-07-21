function [VaR, ETL]=rollingCornishFisher(returns,windowSize, inv0, conf)
% rollingCornishFisher: Function to Calculate Rolling Value-At-Risk
%and Expected Tail Loss Estimates Using a Cornish Fisher Expansion

%This function uses the same logic as all other rolling window estimate
%functions. An estimation window is pulled from the returns series and an
%estimate of VaR and ETL is calculated. Time is then shifted by one day and
%repeated until the final data point in the returns series

%INPUTS:
%       returns: (nx1); Vector of returns.
%       windowSize: (scalar); Size of estimation window for calculating
%                               VaRs
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Investment or Portfolio Size

%OUTPUTS:
%       VaR (nEstimates + 1); Rolling Value-At-Risk estimates at confidence level 'conf'
%       ETL (nEstimates + 1); Rolling Expected Tail Loss estimates at confidence level 'conf'

%***************************************************************************************

nEstimates = length(returns) - windowSize;

VaR = zeros(nEstimates, 1);
ETL =  zeros(nEstimates, 1);

for i = 1:nEstimates + 1
    t = windowSize + i;
%     T = t + 999; 
    
    estWindow = returns(t-windowSize:t-1);
    [VaR(i), ETL(i)] = cornishFisher(estWindow, conf, inv0);
    
end 

