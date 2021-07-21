function CFVAR=CFvar(returns,conf,inv0)
%INPUT
%Rt: returns series
%conf: VaR level
%inv0:  (scalar); Investment or Portfolio Size
%OUTPUT:
%CFVAR: Cornish fisher VaR with alpha probability and time horizon as the returns
%frequancy
% calculate cornish fisher VaR
z=norminv(conf,0,1); %
sigma2=var(returns);
skew=skewness(returns);
kurt=kurtosis(returns);
CFVAR=(mean(returns)+(z+(1/6)*(z^2-1)*skew+(1/24)*(z^3-3*z)*(kurt-3)-(1/36)*(2*z^3-5*z)*skew^2)*sigma2^0.5)*inv0;