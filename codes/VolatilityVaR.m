function [ HWVaR, ETL ] = VolatilityVaR( Returns, inv0, conf, n, sam, lambdaHW)

sigmaHW=zeros(length(n-sam),1);

for t=sam+1:n % Loop to calculate the volatility weighted returns
    
periods=zeros(length(sam),1);
for i=1:sam;
period=((lambdaHW^(i-1)*Returns(t+1-i).^2));  % Calculating the lambda decay factor times squared returns for each sample period in the EWMA
periods(i)=period;
end

Sum=sum(periods);
sigHW=sqrt((1-lambdaHW)*Sum); % Volatility estimate from EWMA for each time period
sigmaHW(t,:)=sigHW;

end

VolWRet=Returns(sam+2:n).*(sigmaHW(end)./sigmaHW(sam+1:n-1)); % Volatility weighted returns based on most recent volatility estimate


HWVaR=historicalVaR(VolWRet, conf, inv0 ) % Volatility Weighted VaR
ETL = historicalETL(VolWRet, conf, inv0 );
end

