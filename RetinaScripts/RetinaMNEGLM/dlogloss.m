%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%%%***###!!!$$$^^^<<<{{{[[[(((%%%)))]]]}}}>>>^^^$$$!!!###***%%%%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%                                                             %%%
%%%               Gradient of the log loss function             %%%
%%%                                                             %%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%%%***###!!!$$$^^^<<<{{{[[[(((%%%)))]]]}}}>>>^^^$$$!!!###***%%%%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%

% Author: Jeff Fitzgerald, August 2010 (updated on Sept 5, 2010)

function df = dlogloss(p, stim, avgs, order)

[Nsamples,Ndim] = size(stim);

ptemp = p(2:Ndim+1);
if order>1
    J = reshape(p(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim]);
end

if order==1
    pSpike = 1./(1+exp(p(1)+stim*ptemp'));  % Nsamples x 1
    averages = mean(pSpike);
    averages(2:Ndim+1,1) = stim'*pSpike/Nsamples;
elseif order==2
    pSpike = 1./(1+exp(p(1)+stim*ptemp'+sum(stim.*(stim*J),2)));  % Nsamples x 1
    averages = mean(pSpike);
    averages(2:Ndim+1,1) = stim'*pSpike./Nsamples;
    temp = stim'*(repmat(pSpike,[1,Ndim]).*stim)./Nsamples;  % Ndim x Ndim
    temp = reshape(temp,[Ndim^2,1]);
    averages(Ndim+2:Ndim+1+Ndim^2) = temp;    
end

df = (avgs - averages)';  % 1 x Ndim