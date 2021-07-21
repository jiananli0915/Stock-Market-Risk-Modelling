%% Read Price Data:
clc;clear
close all
pfe = xlsread('Book1.xlsx');


prices = pfe;
%% Age Weighted VaR:
returns=prices
conf=0.95
inv0=1
% here we do a sensitivity testing
lambda = 0.8; %for age weighting
[ageWeightedVaR99, ageETL99] = ageWeightedEWMA(returns, lambda, inv0, conf);
AgeVar_AMC=AgeweightedHistoricalsimulation(returns,inv0,conf,conf)
%% Historical VaRs Using Tutorial Function:
returns=prices
conf=0.95
inv0=1

[ hVaR99, up2, l2, ETLhist99 ] = historicalVaR( returns, conf, inv0 );

fprintf(['99%% 1-Day Non Parametric VaR:\n<strong>$%f</strong>\n\n'],  hVaR99)
%% Setup:
inv0 = 421328; 
conf = 0.99;
nu = 4;
sigma=0.1417
mu=0.0367

%% Normal 1-Day VaR @ 99%:
%   Returns obviously aren't normal/ i.i.d ~~> Crude estimate
%   Normal VaR is a useful benchmark though.

mu = mean(returns)';
sigma = std(returns);

%formula for Normal VaR in slides
paramVaR99 = -((-mu + sigma*norminv(1-conf,0,1))*inv0);

fprintf(['99%% 1-Day Parametric VaR Under Normal:\n<strong>$%f</strong>\n\n'], paramVaR99)
%% VaR Under Student's t-Distribution With v Degrees of Freedom:
numberOfDegrees = 7; 
tV99 = zeros(numberOfDegrees,1);

for v = 2:numberOfDegrees
    
    tcrit99 = tinv(1 - conf, v);
    tscalar = sqrt((v - 2)/v);
   
    tV99(v-1) = -((-mu + sigma*tscalar*tcrit99)*inv0);
end

%2 DoF fits data best 
tParamVaR99 = tV99(6);

fprintf(['99%% 1-Day Parametric VaR Under T (4 DoF):\n<strong>$%f</strong>\n\n'], tParamVaR99)
