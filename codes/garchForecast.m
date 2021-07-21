function [forecastVec, aic, bic, realizedVol, historicalSigma] = garchForecast(returns, model, nu, P, Q,  trainingSize)
%Function to calculate one-step-ahead volatility forecasts using GARCH:

%Takes a vector of returns and some parameters and returns (k x 1) 
%vectors of GARCH-estimated volatilities, realized volatilities and a
%rolling standard deviation

%Calculates one-step-ahead GARCH forecasts on a rolling basis
%constant sample size of k (see 'INPUTS') -> changing starting date each
%time

%The loop works as follows: 
%   I use a default trainingSize of 500 samples to
%   calibrate a GARCH model for some estimationWindow and use this to calculate a one-step
%   ahead forecast. 


%INPUTS:
%       returns: (nx1); Vector of returns.
%       model: 
%                   MUST BE IN [0,1,2]--> encoded data to specify which GARCH model to use.
%                   0 = GARCH(P,Q)
%                   1 = EGARCH(P,Q)
%                   2 = GJR(P,Q)
%                   3 = GARCH(P,Q) (Under T with nu Degrees of Freedom

%       P: (scalar) Number of autoregressive (ARCH) lags 
%       Q: (scalar); Number of moving average (GARCH) lags 
%       trainingSize: (scalar); Number of observations used to calibrate model
%***************************************************************************************

%BODY:

%Setting default values:
if nargin ==2
    nu = 5;
    P= 1;
    Q = 1; 
    trainingSize = 500;
    
elseif nargin == 3
    P = 1;
    Q = 1; 
    trainingSize = 500;
    
elseif nargin == 4
    Q = 1;
    trainingSize = 500;
    
elseif nargin == 5 
    trainingSize = 500;
      
end 

%Required params:
n = length(returns);
estimates = n - trainingSize;

%Encode GARCH models:
if model == 0
    mdl = garch(P,Q);
elseif model == 1
    mdl = egarch(P,Q);
elseif model == 2
    mdl = gjr(P,Q);
elseif model == 3
    mdl = garch(P,Q);
    mdl.Distribution = struct('Name', 't', 'DoF', nu);
else 
    fprintf('<strong>Wrong Input. Specify Model Using 0 (GARCH), 1 (EGARCH),  2 (GJR) or 3 (T)<strong>')
    return
end
    

%Initializing vectors for recursive loop
forecastVec = zeros(estimates, 1);
historicalSigma = zeros(estimates,1);
aic = zeros(estimates,1);
bic = zeros(estimates,1);


%MAIN LOOP: 
%Calculating Volatility Recursively:

%Starting estimates immediately after initial training data ends:
startingPoint = trainingSize-1;

%Reset Default Algorithm (IP is better for doing loops like this one)
options = optimoptions(@fmincon,'Algorithm','interior-point');


%Calculate Volatility:
for i = 1:estimates + 1
    t = startingPoint + i; %first forecast is the first out-of-training sample
 
    estimationReturns = returns(i:t, 1); %sample of 'trainingSize' to calibrate each model
        
    [est, cov, logL] = estimate(mdl, estimationReturns, 'options', options); %specify model
        
    forecastVec(i) = sqrt(forecast(est, 1, 'Y0', estimationReturns)); %forecast t + 1
    
    %Information Criteria Over Time:
    params = sum(any(cov));   
    [aic(i), bic(i)] = aicbic(logL, params, length(estimationReturns));
    
    %Rolling Standard Deviation:
    historicalSigma(i) = std(estimationReturns);
    
end

%Realized volatility for samples:
realizedVol = sqrt(returns(trainingSize +1:end).^2);

end 

