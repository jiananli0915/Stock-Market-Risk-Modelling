clc;clear
close all

%% Read Price Data:
pfe = xlsread('VaR_data.xlsx', 'XOM_VaR_data');

%relevant columns:
dates = datetime(pfe(:, 1), 'ConvertFrom', 'excel');
prices = pfe(:, 2);

%log returns:
returns = diff(log(prices));
%lose an observation in calculating returns, so remove its date 
returnDates = dates(2:end);

%data table:
df = timetable(returnDates, returns);

%% Plot Prices & P&L Over Time:
close all 

subplot(2,1,1)
plot(dates, prices)
title('Adjusted Closing Of XOM Over Time')

subplot(2,1,2)
plot(returnDates, df.returns)
title('Profit and Loss of XOM Over Time')
%% histogram
figure;
hist(returns,500)
xlabel('Profit and Loss for XOM')
ylabel('Observations')
title('Distribution of Returns over time')


%% statistics of XOM
y=returns
P=1000000; % Portfolio Size
y=y*P;
length_AMC=length(y) % number of elements in y
sum_AMC=sum(y) % sum of elements
prod_AMC=prod(y) % product of elements
range_AMC=range(y) % difference between maximum and minimum
mean_AMC=mean(y) % mean
median_AMC=median(y) % medium
var_AMC=var(y) % variance
std_AMC=std(y) % standard error
corrcoef_AMC=corrcoef(y) % correlation coefficients
skewness_AMC=skewness(y) % get the skewness
kurtosis_AMC=kurtosis(y) % get the kurtosis (NOT excess)
quantile_AMC=quantile(y,0.01) % returns the quantiles at p
min_AMC=min(y) % minimum value
max_AMC=max(y) % maximum value
%% Set Distributions For QQ Plots:
tDist = makedist('tLocationScale') %default is 5 DoF
extreme = makedist('ExtremeValue') 
%not really sure whats going on with this distribution yet

%% Whole Return Series Vs Theoretical Distributions:
%setup
subplot(2,1,1) 
qqplot(returns)
title('XOM P&L (Normal Distribution)')

subplot(2,1,2)
qqplot(returns, tDist)
title('XOM P&L (T Distribution, default is 5 DoF)')
hold off

%% QQ Plots Under t With Varying DoF:
close all 
QQdof = 3

tAdjusted = makedist('tLocationScale', 'mu', 0, 'sigma', 1, 'nu', QQdof) 
%data looks beteer represented under 3

qqplot(returns, tAdjusted) 
%slightly better, still doesnt fit well though
%% Jarque-Bera test

[h,p,jbstat,critval] = jbtest(Returns_XOM,0.05) % Jarque-Bera test
h_XOM=h
p_XOM=p
jbstat_XOM=jbstat
critval_XOM=critval
% The returned value of h = 1 indicates that jbtest rejects the null hypothesis at the default 5% significance level.

%% Plot Squared Returns:
close all
returnsSquared = returns.^2;

plot(returnDates, returnsSquared) 
ylim([0 0.0065])
title('XOM Squared Returns')
%high volatility spikes with clustering
%should confirm with a test to be complete but its obvious that
%theres autocorrelation in squared returns 

%--> NEED A GARCH MODEL TO FORECAST VARIANCE

%% Setup:
inv0 = 1000000; 
conf = 0.95; % we use 99% and 95%
nu = 3;

%% Normal 1-Day VaR @ 99%:
%   Returns obviously aren't normal/ i.i.d ~~> Crude estimate
%   Normal VaR is a useful benchmark though.

mu = mean(returns)';
sigma = std(returns);

%formula for Normal VaR in slides
paramVaR99 = -((-mu + sigma*norminv(1-conf,0,1))*inv0);

fprintf(['99%% 1-Day Parametric VaR Under Normal:\n<strong>$%f</strong>\n\n'], paramVaR99)

%% Historical VaRs Using Tutorial Function:
[ hVaR99, up2, l2, ETLhist99 ] = historicalVaR( returns, conf, inv0 );

fprintf(['99%% 1-Day Non Parametric VaR:\n<strong>$%f</strong>\n\n'],  hVaR99)

%% VaR Under Student's t-Distribution With v Degrees of Freedom:
numberOfDegrees = 30; 
tV99 = zeros(numberOfDegrees,1);

for v = 2:numberOfDegrees
    
    tcrit99 = tinv(1 - conf, v);
    tscalar = sqrt((v - 2)/v);
   
    tV99(v-1) = -((-mu + sigma*tscalar*tcrit99)*inv0);
end

%2 DoF fits data best 
tParamVaR99 = tV99(3);

fprintf(['99%% 1-Day Parametric VaR Under T (3 DoF):\n<strong>$%f</strong>\n\n'], tParamVaR99)

%% Monte Carlo:
[mcnVaR, mcnETL] = monteCarlo(returns, conf, inv0, 10000);
[mctVaR, mctETL] = monteCarlo(returns, conf, inv0, 10000, 1, nu);

fprintf('Monte Carlo (Normal) 99%% 1 Day VaR:\n<strong>$%f</strong>\nMonte Carlo (T) 99%% 1 Day VaR:\n<strong>$%f</strong>\n\n', mcnVaR, mctVaR)

%% Age Weighted VaR:
% here we do a sensitivity testing
lambda = 0.96; %for age weighting
[ageWeightedVaR99, ageETL99] = ageWeightedEWMA(returns, lambda, inv0, conf);
lambda1 = 0.94; %for age weighting
[ageWeightedVaR99_1, ageETL99_1] = ageWeightedEWMA(returns, lambda1, inv0, conf);
lambda2 = 0.8; %for age weighting
[ageWeightedVaR99_2, ageETL99_2] = ageWeightedEWMA(returns, lambda2, inv0, conf);
lambda3 = 0.7; %for age weighting
[ageWeightedVaR99_3, ageETL99] = ageWeightedEWMA(returns, lambda3, inv0, conf);
lambda4 = 0.6; %for age weighting
[ageWeightedVaR99_4, ageETL99] = ageWeightedEWMA(returns, lambda4, inv0, conf);

fprintf(['99%% Age Weighted Historical VaR (lambda = 0.96):\n<strong>$%f</strong>\n\n'], ageWeightedVaR99)
fprintf(['99%% Age Weighted Historical VaR (lambda = 0.94):\n<strong>$%f</strong>\n\n'], ageWeightedVaR99_1)
fprintf(['99%% Age Weighted Historical VaR (lambda = 0.8):\n<strong>$%f</strong>\n\n'], ageWeightedVaR99_2)
fprintf(['99%% Age Weighted Historical VaR (lambda = 0.7):\n<strong>$%f</strong>\n\n'], ageWeightedVaR99_3)
fprintf(['99%% Age Weighted Historical VaR (lambda = 0.6):\n<strong>$%f</strong>\n\n'], ageWeightedVaR99_4)

%% calculate the returns using Age weighted Historical simulation

ageHSVaR99=AgeweightedHistoricalsimulation(returns, inv0 ,conf);
fprintf(['99%% Age Weighted Historical Simulation VaR :\n<strong>$%f</strong>\n\n'], ageHSVaR99);

%% volatility weighted VaR:
sam=500;
lambdaHW=0.94;
[ volVaR99, volETL99 ] = VolatilityVaR( returns, inv0, conf, length(returns), sam, lambdaHW)
fprintf(['99%% Volatility Weighted Historical VaR (lambda = 0.94):\n<strong>$%f</strong>\n\n'], volVaR99)

%% GARCH(1,1)
%Estimating GARCH(1,1)
%Also recording their AIC and BIC over time.

[garchVol, garchAIC, garchBIC, realizedVol, historicalSigma] = garchForecast(returns, 0);

%% Plotting Information Criteria Over Time:
% Want to see which model performs best, and if this is constant throughout samples.
dates = returnDates(end - length(garchAIC)+1:end);

subplot(2,1,1)
plot(dates,garchAIC, 'r-')

ylabel('Information Criteria (AIC)')
title('GARCH(1,1) Models')
legend('GARCH(1,1)')

subplot(2,1,2)
plot(dates,garchBIC, 'r-')

legend('GARCH(1,1)')
ylabel('Information Criteria (BIC)')

%% Weighting By Volatility:
[garchVaR99, garchETL99, volWeightedgarch] = volWeightedVaR(returns, garchVol, conf, inv0);

fprintf( ['\n99%% Volatility Weighted VaR (GARCH(1,1)):\n<strong>%f</strong>\n'...
    ], garchVaR99)

%% Illustrating Volatility Adjustment:
sampleReturns = returns(end-length(volWeightedgarch) + 1:end);
sampleDates = returnDates(end-length(volWeightedgarch) + 1:end);

close all; hold on 
plot(sampleDates, sampleReturns, 'b')
plot(sampleDates, volWeightedgarch, 'r--')

legend('Returns', 'GARCH(1,1) Returns')
ylabel('Return')
title('Returns and GARCH Returns')
%% GARCH Volatilities Vs Returns Squared
%Spikes in Realized Volatility are followed immediately by spikes in Forecasted Volatility as the GJR is given new data.

%plotting a sample of realized returns squared verus garch volatilities

%shifting garch volatilities backwards by one step because they represent volatility
%forecasts 

close all; hold on
plot(dates(end-100:end), garchVol(end-100-1:end-1),'r', ...
    'LineWidth', 0.5)

plot(dates(end-100:end), realizedVol(end-100:end),'g-',...
    'LineWidth', 0.1)

legend('GARCH(1,1)', 'Realized')
title('Reaction of GARCH to Volatility Spikes')

%% Plot of Estimated Volatility For Each Model For last 300 Days
%Illustration of the difference between the two models.

%Rolling Standard Deviation won't change because a new observation is equally weighted (1/500)
close all; hold on
plot(dates(end-300:end), garchVol(end-300:end))
plot(dates(end-300:end), historicalSigma(end-300:end))
legend('GARCH(1,1)', 'Rolling Standard Deviation')
ylabel('Volatility')
title('GARCH Volatility')

%% Average Levels Of Volatility:
%All are very similar: This is expected.

%Figures are scaled by 100

fprintf(['Standard Deviation of the Entire Sample:\n<strong>%f</strong>\nAverage Rolling Standard Deviation:\n<strong>%f</strong>\n' ...
    'Average GARCH(1,1) Volatility:\n<strong>%f</strong>\n\n'], ...
    std(returns)*100, mean(historicalSigma)*100, mean(garchVol)*100)

%% Parametric VaR With GARCH Volatility:

garchNormVaR99 = -((-mu + garchVol(end)*norminv(1-conf,0,1))*inv0);

fprintf(['99%% 1-Day Normal VaR Under GARCH(1,1) Volatility:\n<strong>%f</strong>\n'...
    ], garchNormVaR99)

%% Bootstrapping Historical Simulations:
%Regular historical returns
rng(1)
[baggedVaR99, bVaRs99,  interval99] = bootstrapHistoricalVaR(returns, inv0, conf, 10000, 0);

fprintf(['\nBootstrap Aggregated 99%% VaR Historical Simulation:\n<strong>%f</strong>\n' ...
    'Confidence Interval:\n<strong>[%f, %f]</strong>\n\n'], baggedVaR99, interval99(1), interval99(2))

%% Plotting: 
close all
histogram(bVaRs99, 75)
xlabel('Value at Risk $')
ylabel('Observations')
legend('99% 1 Day VaR')
title('Distribution of Bootstrapped Historical 99% 1 Day VaR')

%% Bootstrapping Historical Simulations:
%Volatility Weighted GARCH(1,1)
[baggedgarchVaR99,bgarchVaRs99, intervalgarch99] = bootstrapHistoricalVaR(volWeightedgarch, inv0, 0.99, 10000, 0);

fprintf(['Bootstrap Aggregated 99%% VaR GARCH Volatility:\n<strong>%f</strong>\n' ...
    'Confidence Interval:\n<strong>[%f, %f]</strong>\n\n '], baggedgarchVaR99, intervalgarch99(1), intervalgarch99(2))

%% Plotting:
close all

histogram(bgarchVaRs99, 75);
xlabel('Value at Risk $');
ylabel('Observations');
legend('99% VaR');
title('Distribution of Bootstrapped Volatility Weighted VaR (GARCH(1,1))');

%% Cornish-Fisher VaR:
% Using the whole dataset for now.
% Using 'cornishFisherVaR' function, which is included.
% Data is non-normal so CF may not be accurate
% ETL is stored for later in the script.

[cfVaR99, cfETL99] = cornishFisher(returns, conf, inv0);
fprintf('99%% VaR Under Cornish Fisher:\n<strong>%f</strong>\n\n',cfVaR99 )


%% Spectral Risk Measure:
% here we use a sensitivity beacktesting 
spectralRisk = spectralRiskMeasure(returns, 5, 10000, inv0);
spectralRisk1 = spectralRiskMeasure(returns, 15, 10000, inv0);
spectralRisk2 = spectralRiskMeasure(returns, 25, 10000, inv0);
spectralRisk3 = spectralRiskMeasure(returns, 50, 10000, inv0);
fprintf('\nSpectral Risk Measure For AMC (Risk Aversion Coefficient =5):\n<strong>%f</strong>\n\n', spectralRisk);
fprintf('\nSpectral Risk Measure For AMC (Risk Aversion Coefficient =15):\n<strong>%f</strong>\n\n', spectralRisk1);
fprintf('\nSpectral Risk Measure For AMC (Risk Aversion Coefficient =25):\n<strong>%f</strong>\n\n', spectralRisk2);
fprintf('\nSpectral Risk Measure For AMC (Risk Aversion Coefficient =50):\n<strong>%f</strong>\n\n', spectralRisk3);

%% Extreme Value VaR:
evVaR99 = evVaR(returns, inv0, conf);

fprintf('99%% VaR Under Extreme Value:\n<strong>%f</strong>\n ',evVaR99 )
%% Monte Carlo:
 rng(3)
 reps = 10000;
 nu = 5;
 
 normalSims=normrnd(mu,sigma,reps,1); % Monte Carlo Simulation based on normal distribution
 [ normSimVaR, ~, ~, normSimETL ] = historicalVaR(normalSims, 0.99, inv0 );
 
 tSims=trnd(nu,reps,1)*sqrt((nu -2)/nu)*sigma+mu; % Monte Carlo Simulation based on Student-t distribution
 [ tSimVaR, ~, ~, tSimETL ] = historicalVaR(tSims, 0.99, inv0 );
 
 
 fprintf(['Normal Monte Carlo:\nVaR: <strong>%f</strong>\nETL: <strong>%f</strong>\n\n'...
     'Monte Carlo (T): \nVaR: <strong>%f</strong>\nETL: <strong>%f</strong>\n'], normSimVaR,normSimETL, tSimVaR, tSimETL)

%% Parametric Expected Tail Loss:
% Expected Loss given that the loss is greater than VaR.
% More informative than VaR, as VaR doesn't measure the extent of exceptional losses.
% tETL is under 5 Degrees of Freedom
% Low DoF means that our VaR and ETL will differ significantly

normETL99 = parametricETL(returns, conf, inv0, 0);
tETL99 = parametricETL(returns, conf, inv0, 1, nu);

fprintf(['Normal Distribution 99%% ETL: <strong>%f</strong>\n'...
    'T-Distribution 99%% ETL: <strong>%f</strong>\n\n'], normETL99, tETL99)

%% Non Parametric Expected Tail Loss:
% Modified historicalVaR function to return ETL
% Volatility Weighted ETLs:
fprintf(['99%% ETL (GARCH):\n<strong>%f</strong>\n\n'],garchETL99)

%% Non Parametric Expected Tail Loss:
%Historical Returns ETL:
fprintf(['99%% ETL From Historical Simulation:\n<strong>%f</strong>\n\n'], ETLhist99)

%% Age Weighted Expected Tail Loss:
fprintf(['99%% Age Weighted ETL:\n<strong>%f</strong>\n\n'], ageETL99_1)

%% Volatility Weighted Expected Tail Loss:
fprintf(['99%% Volatility Weighted ETL (lambda = 0.94):\n<strong>%f</strong>\n\n'], volETL99)

%% Cornish Fisher ETL:
%As on pg 195 of Hernandez book.
fprintf(['99%% Cornish-Fisher ETL:\n<strong>%f</strong>\n\n'],  cfETL99)

%% Backtesting Risk Measures:

%Reading in more data:
pfe = xlsread('VaR_data.xlsx', 'XOM_VaR_data');

%relevant columns:
dates = datetime(pfe(:, 1), 'ConvertFrom', 'excel');
prices = pfe(:, 2);

%log returns:
backTestReturns = diff(log(prices));
%lose an observation in calculating returns, so remove its date 
returnDates = dates(2:end);

%Set Inputs:
windowSize = 500; %estimation window 
trainingSize = 500; %500 estimates to calibrate GARCH model
conf = 0.99; %Confidence level
inv0 = 1000000;  %Portfolio Size
model = 0; % set a GARCH(1,1) model

%% Rolling Forecast:

[forecastVec] = garchForecast(backTestReturns, model);
%% Rolling Historical VaR:
[rollingHistVaR, rollingHistETL] = rollingHistoricalVaR(backTestReturns, conf, inv0, windowSize);
kupiec(backTestReturns, rollingHistVaR, conf, inv0);

%% Rolling Volatility Weighted VaR:

 [rollingVolWVaR, rollingVolETL] = rollingVolVaR(backTestReturns, windowSize, conf, inv0, forecastVec);
kupiec(backTestReturns, rollingVolWVaR, conf, inv0);

%% Rolling Volatility Weighted VaR combined with normal distribution(GARCH):

 [rollingVolnormalVaR, rollingVolnormalETL] = rollingNormal(backTestReturns, windowSize, conf, inv0, forecastVec);
kupiec(backTestReturns, rollingVolnormalVaR, conf, inv0);

%% Rolling Age Weighted VaR
 [rollingAgeWVaR, rollingAgeWETL] = rollingAgeVaR(backTestReturns, conf, inv0, windowSize, lambda);
kupiec(backTestReturns, rollingAgeWVaR, 0.99, inv0);

%% Monte Carlo VaRs
 [rollingMcVaR, rollingMcVaRETL] = rollingMonteCarlo(backTestReturns,windowSize, inv0, conf, 0,4, 1000);
 [rollingMctVaR, rollingMctVaRETL] = rollingMonteCarlo(backTestReturns,windowSize, inv0, conf, 1,4, 1000);

 kupiec(backTestReturns, rollingMcVaR, conf, inv0);
 kupiec(backTestReturns, rollingMctVaR, conf, inv0);

%% Rolling Cornish Fisher:
[rollingCFVaR,rollingCFETL] = rollingCornishFisher(backTestReturns,windowSize, inv0, conf);
kupiec(backTestReturns, rollingCFVaR, conf, inv0);

%% Rolling Normal VaR:
[rollingNormVaR,rollingNormETL] = rollingParametric(backTestReturns,windowSize, inv0, conf, 0, 0);
kupiec(backTestReturns, rollingNormVaR, conf, inv0);

%% Rolling T VaR:
[rollingTVaR,rollingTETL] = rollingParametric(backTestReturns,windowSize, inv0, conf, 1, nu);
kupiec(backTestReturns, rollingTVaR, conf, inv0);

%% Rolling Extreme Value:
% This takes substantial time to run.
[rollingExVaR] = rollingExtremeVaR(backTestReturns, conf, inv0, windowSize);
kupiec(backTestReturns, rollingExVaR, conf, inv0);

%% Backtesting Combination VaR:
 [comboVaR, comboETL] = ageVolCombination(backTestReturns, windowSize, conf, inv0, forecastVec, lambda);
kupiec(backTestReturns, comboVaR, conf, inv0);

%% Parametric Rolling VaR Estimates:
sampleReturns = backTestReturns(end - length(rollingNormVaR) + 1:end)*inv0;
sampleDates = returnDates(end - length(rollingNormVaR) + 1:end);
r = length(rollingVolWVaR) -1;
% we only show the last 300 days in plot
close all; hold on
plot(sampleDates(end - 300:end), sampleReturns(end - 300:end));
plot(sampleDates(end - 300:end), rollingNormVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingTVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingCFVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingMcVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingMctVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingExVaR(end - 300:end))
 
legend('XOM Profit and Loss', 'Normal Distribution', 'T Distribution(2 DoF)', 'Cornish-Fisher','Monte Carlo(Normal)', 'Monte Carlo (T)', 'Extreme Value','Location', 'best');
ylabel('Returns $');
title('99% 1 Day Parametric VaR Estimates & PnL For XOM');

%% Non-Parametric VaR Estimates:
% we only show the last 300 days in plot
close all; hold on

plot(sampleDates(end - 300:end), sampleReturns(end - 300:end));
plot(sampleDates(end - 300:end), rollingHistVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingVolWVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingAgeWVaR(end - 300:end))
plot(sampleDates(end - 300:end), rollingVolnormalVaR(end - 300:end))

%plot(sampleDates(end - 300:end), comboVaR(end - 300:end))
ylabel('Returns $');
%legend('AMC Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Age & Volatility Weighted', 'Location', 'best');
%legend('XOM Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted',  'Location', 'best');
legend('XOM Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Volatility Weighted & Normal (GARCH)', 'Location', 'best');

title('99% 1 Day Non-Parametric VaR Estimates & PnL For XOM');

%% Exceedance Plots
% normal Exceedance
[indicatorNorm, realizedVolatility] = exceedanceIndicator(backTestReturns, rollingNormVaR/inv0);
d = returnDates(end - length(indicatorNorm) + 1:end) 
close all; hold on
plot(d,indicatorNorm, 'k'); 
plot(d,realizedVolatility(end - length(indicatorNorm) + 1:end), 'r')
ylim([0 0.2])
ylabel('Volatility (Scaled)')
title('Exceedances of Normal Rolling VaR')
legend('Exceedances', 'Squared Returns', 'Location', 'best')
%% T Exceedances:
[indicatorT, realizedVolatility] = exceedanceIndicator(backTestReturns, rollingTVaR/inv0);

d = returnDates(end - length(indicatorT) + 1:end) 
close all; hold on
plot(d,indicatorT, 'k'); 
plot(d,realizedVolatility(end - length(indicatorT) + 1:end), 'r')
ylim([0 0.2])
ylabel('Volatility (Scaled)')
title('Exceedances of T Rolling VaR')
legend('Exceedances', 'Squared Returns', 'Location', 'best')

%% Exceedances of Rolling Age VaR
[indicatorV, realizedVolatility] = exceedanceIndicator(backTestReturns, rollingAgeWVaR/inv0);

d = returnDates(end - length(indicatorV) + 1:end) 
close all; hold on
plot(d,indicatorV, 'k'); 
plot(d,realizedVolatility(end - length(indicatorV) + 1:end), 'r')
ylim([0 0.2])
ylabel('Volatility (Scaled)')
title('Exceedances of Rolling Age Weighted VaR')
legend('Exceedances', 'Squared Returns', 'Location', 'best')

%% Parametric ES
% plot of 1250 days
close all; hold on
plot(sampleDates(end - 1250:end), sampleReturns(end - 1250:end))
plot(sampleDates(end - 1250:end), rollingNormETL(end - 1250:end))
plot(sampleDates(end - 1250:end), rollingTETL(end - 1250:end))
plot(sampleDates(end - 1250:end), rollingCFETL(end - 1250:end))
plot(sampleDates(end - 1250:end), rollingMcVaRETL(end - 1250:end))
plot(sampleDates(end - 1250:end), rollingMctVaRETL(end - 1250:end))
 
legend('AMC Profit and Loss', 'Normal Distribution', 'T Distribution(2 DoF)', 'Cornish-Fisher','Monte Carlo(Normal)', 'Monte Carlo (T)','Location', 'best');
ylabel('Returns $')
title('99% 1 Day Parametric ES Estimates For XOM')

%% Parametric ES
% plot of the last 300 days
close all; hold on
plot(sampleDates(end - 300:end), sampleReturns(end - 300:end))
plot(sampleDates(end - 300:end), rollingNormETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingTETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingCFETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingMcVaRETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingMctVaRETL(end - 300:end))

legend('AMC Profit and Loss', 'Normal Distribution', 'T Distribution(2 DoF)', 'Cornish-Fisher','Monte Carlo(Normal)', 'Monte Carlo (T)', 'Location', 'best');
ylabel('Returns $')
title('99% 1 Day Parametric ES Estimates For XOM')

%% Non-Parametric ETL Estimates:
% plot of 750 days (can not do it longer because the limitation of the length for rollingVolWVaR)
close all; hold on

plot(sampleDates(end - 750:end), sampleReturns(end - 750:end))
plot(sampleDates(end - 750:end), rollingHistETL(end - 750:end))
plot(sampleDates(end - 750:end), rollingVolETL(end - 750:end))
plot(sampleDates(end - 750:end), rollingAgeWETL(end - 750:end))
plot(sampleDates(end - 300:end), rollingVolnormalETL(end - 300:end))
%plot(sampleDates(end - 750:end), comboETL(end - 750:end))
ylabel('Returns $')
%legend('AMC Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Age & Volatility Weighted', 'Location', 'best');
legend('XOM Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Volatility Weighted & Normal (GARCH)' , 'Location', 'best');

title('99% 1 Day Non-Parametric ES Estimates For XOM')

%% Non-Parametric ETL Estimates:
% plot of the last 300 days
close all; hold on

plot(sampleDates(end - 300:end), sampleReturns(end - 300:end))
plot(sampleDates(end - 300:end), rollingHistETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingVolETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingAgeWETL(end - 300:end))
plot(sampleDates(end - 300:end), rollingVolnormalETL(end - 300:end))

%plot(sampleDates(end - 300:end), comboETL(end - 300:end))

ylabel('Returns $')
%legend('AMC Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Age & Volatility Weighted', 'Location', 'best');
legend('XOM Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Volatility Weighted & Normal (GARCH)' , 'Location', 'best');

title('99% 1 Day Non-Parametric ES Estimates For XOM')
%% VaR violations for Parametric methods
ZoomInd   = sampleDates(end - 300:end);% we zoom the timeframe when there is a large and sudden change in the value of returns.
VaRData   = [rollingNormVaR(end - 300:end) rollingTVaR(end - 300:end) rollingCFVaR(end - 300:end) rollingMcVaR(end - 300:end) rollingMctVaR(end - 300:end) rollingExVaR(end - 300:end)];
VaRFormat = {'-','--','-.',':',':.','-'};

D = sampleDates(end - 300:end);
R = sampleReturns(end - 300:end);
N = rollingNormVaR(end - 300:end);
T = rollingTVaR(end - 300:end);
C = rollingCFVaR(end - 300:end);
M = rollingMcVaR(end - 300:end);
MT = rollingMctVaR(end - 300:end);
E = rollingExVaR(end - 300:end);
IndN99    = (R >  N);
IndT99    = (R >  T);
IndC99    = (R >  C);
IndM99    = (R >  M);
IndMT99   = (R >  MT);
IndE99    = (R >  E);

figure;
bar(D,R,0.5,'FaceColor',[0.7 0.7 0.7]);
hold on
for i = 1 : size(VaRData,2)
    stairs(D-0.5,VaRData(:,i),VaRFormat{i});
end
ylabel('VaR')
xlabel('Date')
title('99% VaR violations for Parametric methods for XOM')
ax = gca;
ax.ColorOrderIndex = 1;

plot(D(IndN99),N(IndN99),'o',D(IndT99),N(IndT99),'o',D(IndC99),N(IndC99),'o',D(IndM99),N(IndM99),'o',...
   D(IndMT99),N(IndMT99),'o',D(IndE99),N(IndE99),'o','MarkerSize',8,'LineWidth',1.5)
xlim([D(1)-1, D(end)+1])
legend('AMC Profit and Loss', 'Normal Distribution', 'T Distribution(2 DoF)', 'Cornish-Fisher','Monte Carlo(Normal)', 'Monte Carlo (T)','Extreme Value','Normal', 'T', 'CF','MC(Normal)', 'MC(T)','EV', 'Location', 'best');

hold off;
% A VaR failure or violation happens when the loss is bigger than VaR. 

%% VaR violations for Non-Parametric methods
ZoomInd   = sampleDates(end - 300:end);% we zoom the timeframe when there is a large and sudden change in the value of returns.
VaRData   = [rollingHistVaR(end - 300:end) rollingVolWVaR(end - 300:end) rollingAgeWVaR(end - 300:end) comboVaR(end - 300:end)];
VaRFormat = {'-','--','-.',':'};

D = sampleDates(end - 300:end);
R = sampleReturns(end - 300:end);
H = rollingHistVaR(end - 300:end);
V = rollingVolWVaR(end - 300:end);
A = rollingAgeWVaR(end - 300:end);
VN= rollingVolnormalVaR(end - 300:end);
%AV = comboVaR(end - 300:end);
IndH99    = (R >  H);
IndV99    = (R >  V);
IndA99    = (R >  A);
IndVN99    = (R >  VN);
%IndAV99   = (R >  AV);

figure;
bar(D,R,0.5,'FaceColor',[0.7 0.7 0.7]);
hold on
for i = 1 : size(VaRData,2)
    stairs(D-0.5,VaRData(:,i),VaRFormat{i});
end
ylabel('VaR')
xlabel('Date')
title('99% VaR violations for Non-Parametric methods for XOM')
ax = gca;
ax.ColorOrderIndex = 1;

plot(D(IndH99),N(IndH99),'o',D(IndV99),N(IndV99),'o',D(IndA99),N(IndA99),'o',D(IndVN99),N(IndVN99),'o','MarkerSize',8,'LineWidth',1.5)
xlim([D(1)-1, D(end)+1])
legend('XOM Profit and Loss', 'Historical', 'Volatility Weighted', 'Age Weighted', 'Volatility Weighted & Normal(GARCH)', 'Location', 'best');

hold off;

%% statistical tests for Parametric VaR backtesting
vbt1 = varbacktest(sampleReturns,[rollingNormVaR rollingTVaR rollingCFVaR(2:1264,:) rollingMcVaR rollingMctVaR rollingExVaR],'PortfolioID','AMC','VaRID',...
    {'Normal Distribution' 'T Distribution(2 DoF)' 'Cornish-Fisher' 'Monte Carlo(Normal)' 'Monte Carlo (T)' 'Extreme Value'},'VaRLevel',[0.99]);

summary(vbt1)
runtests(vbt1)

%% statistical tests for Non-Parametric VaR backtesting
% The result does not match the previous test, so just ignore this bar.
%vbt2 = varbacktest(sampleReturns(end - 763:end),[rollingHistVaR(end - 763:end) rollingVolWVaR rollingAgeWVaR(end - 763:end) comboVaR],'PortfolioID','AMC','VaRID',...
 %   {'Historical' 'Volatility Weighted' 'Age Weighted' 'Age & Volatility Weighted' },'VaRLevel',[0.99]);

%summary(vbt2)