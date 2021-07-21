function [modelCorrect,failureFreq,freqInterval ] = kupiec(returns, VaR, conf, inv0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%make sure that VaR and Returns are the same length:
VaR_matched = VaR(1:end-1);
n = length(VaR_matched);
 
sampleReturns = returns(end - n + 1:end);

%get profit & loss:
PnL = inv0*sampleReturns;

% Probability of an exceedance
p=1-conf;

%losses:
tail_losses=-PnL(-PnL>VaR_matched);

% Number of exceedances
x=length(tail_losses) ;    

%test:
modelCorrect=1-binocdf(x,n,p);
[failureFreq,freqInterval]=binofit(x,n,p);

if modelCorrect >= p
    fprintf('\n\nExpected Exceedances: <strong>%d</strong>\nModel Exceedances: <strong>%d</strong>\nViolation Ratio: <strong>%f</strong>\n\nConclusion: <strong>Fail to Reject The Null Hypothesis (p(model) = %f)</strong>; Model is Acceptable\n\n', round(n*p), x, (x/(round(n*p))), modelCorrect) 
elseif modelCorrect < p
    fprintf('\n\nExpected Exceedances: <strong>%d</strong>\nModel Exceedances: <strong>%d</strong>\nViolation Ratio: <strong>%f</strong>\n\nConclusion: <strong>Reject The Null Hypothesis (p(model) = %f); Model is Rejected</strong>\n\n', round(n*p), x, (x/(round(n*p))), modelCorrect)
end 

end 

