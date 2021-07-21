function [VaR] = evVaR(returns,inv0, conf)
%Function to Calculate VaR Using Extreme Value Theory

%Based on Tutorial Code.

%INPUTS:
%       returns: (nx1); Vector of Returns.
%       inv0: (scalar); Investment or Portfolio Size
%       conf: (scalar); Confidence Level for VaR (i.e 0.95, 0.99)

%OUTPUTS:
%       VaR: (scalar); VaR Estimate

%***************************************************************************************

k=1;
n=25;
steps=20;

MinimalRet=zeros(steps,1);
for i=1:steps:(n*steps)
    Trets=returns(i:i+steps-1);
    MinimalRet(k,:)=min(Trets);
    k=k+1;
end

parmhat = gevfit(MinimalRet);
VaR=(-parmhat(3)*ones(size(conf))+(parmhat(2)/parmhat(1))*(1-(-log(conf^n)).^parmhat(1)))*inv0;      

end 

