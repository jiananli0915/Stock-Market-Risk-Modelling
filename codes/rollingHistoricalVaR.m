function [VaR, ETL] = rollingHistoricalVaR(returns, conf, inv0, estWindowSize)
%rollingHistoricalVaR: Function to calculate rolling historical VaR
    %Takes returns, confindence level, portfolio size and (option) number
    %of samples to use in Estimation Window and returns vector of rolling
    %historical Value-At-Risk estimates
    %Uses 'historicalVaR' function
 
%INPUTS:
%       returns: (nx1): Vector of Portfolio profit & loss
%       conf: (scalar); Confidence Level for VaR & ETL (i.e 0.95, 0.99)
%       inv0: (scalar); Investment or Portfolio Size

%OUTPUTS:
%       VaR (nEstimates x 1): Rolling Value-At-Risk estimates
%******************************************************************************


%Set default value for samples in Estimation Window:
if nargin == 3
    estWindowSize = 1000;
elseif nargin == 4
    estWindowSize = estWindowSize;
end

%Define nEstimates:
nEstimates = length(returns) - estWindowSize;

%Initialize VaR vector:
VaR = zeros(nEstimates, 1);
ETL = zeros(nEstimates, 1);

%Loop to calculate rolling historical VaR
for i = 1:nEstimates + 1
    t = estWindowSize + i - 1; 
    estWindow = returns(i:t);
    VaR(i) = historicalVaR(estWindow, conf, inv0);
    ETL(i) =  historicalETL(estWindow, conf, inv0);
end 



