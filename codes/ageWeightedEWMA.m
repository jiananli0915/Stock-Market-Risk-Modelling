function [VaR, ETL, sortStamped] = ageWeightedEWMA(returns, lambda, inv0, conf)
%Function to calculate an Age-Weighted VaR & ETL Using the EWMA Model

%   Takes an nx1 vector of returns, a decay factor & a threshold and
%   returns an Age-Weighted VaR and ETL

%INPUTS:
%       returns: (nx1); Vector of returns.
%       lambda: (scalar); Decay Factor.
%       thres: (scalar); Threshold for VaR (0.01 or 0.05)
%       inv0: (scalar); Investment (portfolio size)
     
%OUTPUTS:
%       ageWeightedVaR: (scalar); Age Weighted VaR Estimate
%       ETL: (scalar); Age Weighted ETL Estimate

%***************************************************************************************

%Starting Params:   
n = length(returns);
flipped = flipud(returns); %most recent at the top
thres = 1-conf; %alpha is our threshold

%Initialize Vector:
lambdaWeight = zeros(n,1);

%Exponentially Weight:
 for i = 1:n
     lambdaWeight(i,1) = ((lambda^(i-1)) * (1 - lambda))/(1 - lambda^n); %EWMA weighting 
 end 

%matrix of returns stamped with their weight (exponentially decreasing):
stamped = [flipped lambdaWeight]; 

%sort matrix by column 1 (smallest , i.e most negative returns are first):
sortStamped = sortrows(stamped,1); 

%initialize cumulative weight variable:
cumulativeWeight = sortStamped(1,2); 

%Main loop to calculate VaR estimate
for j = 2:n
   
    %Cumulative sum of weights:
    cumulativeWeight(j) = cumulativeWeight(j-1) + sortStamped(j,2); 
    
    %If we have a sudden large loss that, when stamped by weight, exceeds
    %our threshold, this is the VaR estimate.
    if cumulativeWeight(1) > thres
        VaR = -sortStamped(1,1)*inv0;

        %Since we don't have mean losses in excess of VaR, set ETL to zero
        %and handle in rolling function:
        ETL = 0;
        break
        
     %if our cumulative weights add up to exactly the threshold value, 
     %the loss at this index is the VaR estimate:
    elseif  cumulativeWeight(j) == thres
        idx = j;
        VaR = -sortStamped(idx , 1)*inv0;
        ETL = mean(sortStamped(end-idx+1:end, 1))*inv0;
        break 
    
      %if the cumulative index breaches the threshold, interpolate between 
      %the losses corresponding to this and the previous index.
      %This interpolated value is the VaR estimate.
    elseif cumulativeWeight(j) > thres 
        
        %index & value of lower VaR: 
        lowerIndex = j - 1;
        lowerVaR = sortStamped(lowerIndex, 1);
        lowerWeight = cumulativeWeight(lowerIndex);
        
        %index & value of upper Var: 
        upperIndex = j;
        upperWeight = cumulativeWeight(upperIndex);
        upperVaR = sortStamped(upperIndex,1);
        
        %Get VaR using linear interpolation:
        numerator = lowerVaR*(thres - lowerWeight)  + upperVaR*(upperWeight - thres);
        denominator = upperWeight - lowerWeight;

        VaR= -(numerator/denominator)*inv0;
        ETL = mean(sortStamped(end-upperIndex+1:end, 1))*inv0;
        break
    end 
end 

end 



