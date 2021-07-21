function [VaR]=AgeweightedHistoricalsimulation(Return,P,cl,doc)
% P=investment size
% cl=Chosen confidence level
% doc= Degree of confidence
PL=Return*P; % Historical P&L for investment at that size
mu=mean(Return); % Average returns of investment
sigma=std(Return); % Standard Deviation for the returns of investment

%% Set up the loss data as we need
LossesNegative=-PL;   % Transforms losses to positive values for investment
n=length(LossesNegative);

%% Set up lossmatrix and cumulative age-weights
num= 1:n;
num=num'; % Nummber for each scene
lambda=1-(1-cl)/2; % lambda value for this situation
age_weight=zeros(n,1); % age_weight value for each scene
for i=1:n
    age_weight(i) = (lambda^(n-i)*(1-lambda))/(1-lambda^n);
end

Lossesmatrix=[num,LossesNegative,age_weight]; % lossmatrix
Lossesmatrix=sortrows(Lossesmatrix,-2); % lossmatrix
age_weight=Lossesmatrix(:,3);
cum_age_weight=zeros(n,1); % Original cumulative age-weights
for i=2:n % Cumulative age-weights vector
    cum_age_weight(1) = age_weight(1);
    cum_age_weight(i) = age_weight(i)+cum_age_weight(i-1);
end
%% Find out age-weightes VaR
No_ageweight=find(cum_age_weight>=(1-cl),1,'first'); % The first cumulative age-weight higher than alpha
VaR=Lossesmatrix(No_ageweight,2); % VaR value

end