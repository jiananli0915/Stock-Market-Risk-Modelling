function [SRM] = spectralRiskMeasure(returns, R, n, inv0)
%Function to Calculate Spectral Risk Measure of an Investment
%Based on Tutorial Code

%INPUTS:
%       returns: (nx1); Vector of returns.
%       R: (scalar); Risk Aversion Coefficient 
%       n: (scalar); Integration Limit 
%       inv0: (scalar); Investment or Portfolio Size

%OUTPUTS:
%       SRM: (scalar); Spectral Risk Measure

%***************************************************************************************
%BODY:


PL=returns*inv0; % P&L for portfolio

mu=mean(returns); % Average returns of portfolio
sigma=std(returns); % Standard Deviation for the returns of portfolio

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

phi=zeros(n,1); 
var=zeros(n,1); 
f=zeros(n,1); 
for i=1:n
    phi(i)= R*exp(-R*(1-p(i)))/(1-exp(-R)); 
    var(i)=mu+sigma*norminv(p(i),0,1); 
    f(i)=var(i)*phi(i); 
end

%% Estimate the value of SRM
TrapWeightedVals=w'*f;
SRM=TrapWeightedVals*inv0;
end 