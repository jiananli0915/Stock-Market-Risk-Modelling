function [baggedVaR, bootVaR, confInterval] = bootstrapHistoricalVaR(returns, inv0, conf, nboots, type)
%Function to Calculate Bootstrap-Aggregated Historical VaR With Confidence Interval

%INPUTS:
%       returns: (nx1); Vector of Returns.
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Confidence Level for VaR (i.e 0.95, 0.99)
%       nboots: (scalar); Number of Bootstrap Repetitions

%OUTPUTS:
%       bootVaR: (nboots x 1); Vector of VaR Estimates Over nboots samples
%       baggedVaR: (scalar); Bootstrap-Aggregated VaR Estimate
%       confInterval: (1x2); Confidence Interval for VaR Estimates.

%***************************************************************************************


%set seed for reproducibility
rng(1)

%Set up Function Handle to Call 'historicalVaR' Function
if type == 0
    fun = @(returns, inv0, conf ) historicalVaR(returns, conf, inv0); 
elseif type == 1
    fun = @(returns, inv0, conf ) historicalETL(returns, conf, inv0); 
end 
    
%Bootstrap Estimates:
bootVaR = bootstrp(nboots, fun, returns, inv0, conf);

%Average VaR is our Estimate:
baggedVaR = mean(bootVaR);
    
%Confidence Interval:
[confInterval] = bootci(nboots, fun, returns, inv0, conf);

end

