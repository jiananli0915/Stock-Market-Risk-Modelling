function [VaR, ETL] = monteCarlo(returns, conf, inv0, reps, dist, nu)
%monteCarlo: Function To Calculate VaR and ETL Using A Crude Monte Carlo 

%INPUTS:
%       returns: (nx1); Vector of returns.
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Confidence Level for VaR (i.e 0.95, 0.99)
%       reps: (scalar); No. of Replications
%       dist: 0 if Normal, 1 if T
%       nu: Degrees of Freedom
rng(2)
if nargin == 4
    dist = 0; 
end 

mu = mean(returns);
sigma = std(returns);

if dist == 0
    normalMCReturns=normrnd(mu,sigma,reps,1); % Monte Carlo Simulation based on normal distribution
    
    VaR=historicalVaR(normalMCReturns, conf, inv0);
    ETL = historicalETL(normalMCReturns, conf, inv0);
    
elseif dist == 1
    Veps = trnd(nu, reps, 1);
    tscalar = sqrt((nu-2)/2);
    tMCReturns=Veps*tscalar*sigma+mu; % Monte Carlo Simulation based on Student-t distribution
    
    VaR=historicalVaR(tMCReturns, conf, inv0);
    ETL = historicalETL(tMCReturns, conf, inv0);
    
end 
end

    



