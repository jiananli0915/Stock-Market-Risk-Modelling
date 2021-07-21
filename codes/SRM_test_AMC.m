clc;clear
close all

%% Read Price Data:
pfe = xlsread('VaR_data.xlsx', 'AMC_VaR_data');

%relevant columns:
dates = datetime(pfe(:, 1), 'ConvertFrom', 'excel');
prices = pfe(:, 2);

%log returns:
returns = diff(log(prices));
%lose an observation in calculating returns, so remove its date 
returnDates = dates(2:end);

%data table:
df = timetable(returnDates, returns);

%% Setup:
inv0 = 1000000; 
conf = 0.99;
nu = 4;
R1=5;
R2=15;
R3=25;
R4=50;
n=10000;

%% Trapezodial Rule

a=1/n; 
b=(n-1)/n;  
h=(b-a)/(n-1);   

p=zeros(n,1); 
for i=1:n          
    p(i)=a+(i-1)*h; 
end

w=zeros(n,1); 
w(1)=h/2; 
w(n)=h/2;   
for i=2:n-1  
    w(i)=h;
end

%% Specify f(x)for the Spectral Risk Measure 
mu = mean(returns)';
sigma = std(returns);
phi=zeros(n,1); 
var=zeros(n,1); 
f=zeros(n,1); 
for i=1:n
    phi1(i)= R1*exp(-R1*(1-p(i)))/(1-exp(-R1)); 
    var(i)=mu+sigma*norminv(p(i),0,1); 
    f1(i)=var(i)*phi1(i); 
    phi2(i)= R2*exp(-R2*(1-p(i)))/(1-exp(-R2)); 
    f2(i)=var(i)*phi2(i); 
    phi3(i)= R3*exp(-R3*(1-p(i)))/(1-exp(-R3)); 
    f3(i)=var(i)*phi3(i); 
    phi4(i)= R4*exp(-R4*(1-p(i)))/(1-exp(-R4)); 
    f4(i)=var(i)*phi4(i); 
end

%% Estimate the value of SRM
TrapWeightedVals=w'*f;
SRM=TrapWeightedVals*inv0;
%% 
figure;
h1=cdfplot(f1*inv0);
set(h1, 'LineWidth',2); 
hold on
h2=cdfplot(f2*inv0)
set(h2, 'LineWidth',2); 
hold on
h3=cdfplot(f3*inv0)
set(h3, 'LineWidth',2); 
hold on
h4=cdfplot(f4*inv0)
set(h4, 'LineWidth',2); 
xlabel('value of loss function')
ylabel('Cumulative Probability')
legend({'R = 5','R=15','R=25','R=50'})
title('Spectral Risk Measure for AMC')
hold off
