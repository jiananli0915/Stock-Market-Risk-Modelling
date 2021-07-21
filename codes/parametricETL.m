function ETL = parametricETL(returns, conf, inv0, dist, nu)
%Function to Calculate Expected Tail Loss For Normal and T Distributions


%INPUTS:
%       returns: (nx1); Vector of Returns
%       conf: (scalar); Confidence Level for ETL (i.e 0.95, 0.99)
%       inv0: (scalar); Investment of Portfolio Size
%       dist: (scalar); Distribution of ETL
%                            0: Normal Distribution
%                            1: t Distribution
%       nu: (scalar); Degrees of Freedom for T (default is 5)

%OUTPUTS:
%       ETL: (scalar); Expected Tail Loss Under Specified Distribution and
%                           Confidence Level

%***************************************************************************************

%BODY:


%Mean and Standard Deviation of Returns:
mu=mean(returns);
sigma=std(returns);

%Quantile:
a = 1 - conf;

%ETL for Specified Distributions:
%1) Normal
if dist == 0 
    ETLa = (1/a)*normpdf(norminv(conf,0,1), 0, 1);

    ETL = (ETLa*sigma - mu)*inv0;
   
 %2) T
elseif dist == 1 
    xanu = tinv(a, nu);

    ETL = ((1/a)*(1/(nu - 1))* (nu - 2 + xanu^2)*tpdf(xanu, nu)*sigma - mu)*inv0;

end 

end 


