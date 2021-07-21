function [VaR,ETL] = rollingMonteCarlo(returns,windowSize, inv0, conf, dist, nu, reps)
%rollingMonteCarlo: Function to Calculate Rolling Value-At-Risk and
%Expected Tail Loss Estimates Using a Crude Monte Carlo

%INPUTS:
%       returns: (nx1); Vector of returns.
%       windowSize: (scalar); Size of estimation window for calculating
%                               VaRs
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Investment or Portfolio Size
%       dist: Code for Distribution under Which to Estimate VaR & ETL
%                    0 == Normal
%                    1 == T
%       nu: Degrees of Freedom if Estimating Under t (Default is 5)
%       reps: Number of Replications

%OUTPUTS:
%       VaR (nEstimates x 1); Rolling Value-At-Risk estimates at confidence level 'conf'
%       ETL (nEstimates x 1); Rolling Expected Tail Loss estimates at confidence level 'conf'

%Default degrees of Freedom:
if nargin == 5
    nu = 5;
end 

rng(2)
%Number of Estimates:
nEstimates = length(returns) - windowSize;

%Initialize Vectors:
VaR = zeros(nEstimates, 1);
ETL = zeros(nEstimates, 1);

%Normal Distribution:
if dist == 0
    for i = 1:nEstimates 
        
     %mapping for returns vector
        t = windowSize + i;
        
      %pull a sample of size 'windowSize' from returns vector
        estWindow = returns(t-windowSize + 1:t);
    
   
      %Rolling VaR & ETL estimates
        [VaR(i), ETL(i)] =  monteCarlo(estWindow, conf, inv0, reps);

    end 

elseif dist == 1
    for i = 1:nEstimates 
      %mapping for returns vector 
        t = windowSize + i;
        
      %pull a sample of size 'windowSize' from returns vector
        estWindow = returns(t-windowSize+1:t);
        
        [VaR(i), ETL(i)] =  monteCarlo(estWindow, conf, inv0, reps, 1, nu);
        
    
    end 
end 