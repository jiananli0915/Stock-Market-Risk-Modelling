function [ VaR, upper_var, lower_var, ETL ] = historicalVaR(returns, conf, inv0 )
%Calculates Historical VaR based on Tutorial Code
%Modification for ETL:
%           - ETL is defined as the mean of the losses in excess of VaR
%           index
%           - it does not include the VaR index
%           - if we are interpolating, we take the index after upper_index
%           - will be double counting an observation otherwise
%           - cant take a fraction of upper_index either -> will drag mean
%           down


PL    = returns*inv0;                         
LossesNegative = -PL;                       
LossesSorted   = sort(LossesNegative);      
n = length(LossesSorted);                      

index = conf*n  ;                     

%if index is an integer
%ETL is the average of losses in excess of VaR
if index-round( index ) == 0               
   VaR = LossesSorted( index ); 
   upper_var = VaR; lower_var = VaR;
   ETL = mean(LossesSorted(index + 1:end));

% handle non integer index:
%take next value from upper_index for ETL
elseif index-round( index ) > 0 || index-round( index ) < 0 

    upper_index = ceil( index );                
    upper_var   = LossesSorted( upper_index );  

    lower_index = floor( index );              
    lower_var   = LossesSorted( lower_index ); 
    
    ETL = mean(LossesSorted(upper_index + 1 :end));
% If lower and upper indices are the same, VaR is upper VaR
        if upper_index == lower_index
            VaR = upper_var;
        
        else lower_weight = ( upper_index-index )/( upper_index-lower_index );   % weight on lower_var
            upper_weight = ( index-lower_index )/( upper_index-lower_index );   % weight on upper_var
            VaR = lower_weight*lower_var + upper_weight*upper_var;              % VaR
        end    
end



end

