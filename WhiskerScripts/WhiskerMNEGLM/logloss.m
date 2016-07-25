%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%%%***###!!!$$$^^^<<<{{{[[[(((%%%)))]]]}}}>>>^^^$$$!!!###***%%%%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%                                                             %%%
%%%                     The log loss function                   %%%
%%%                                                             %%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%
%%%%%***###!!!$$$^^^<<<{{{[[[(((%%%)))]]]}}}>>>^^^$$$!!!###***%%%%%
%%%%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%   %%%%%

% The log loss is defined as
%  L = -<Log[P(y|x)]>, averaged over the empirical distribution
%  L = -1/Nsamples * sum_Stim(NspikesperStim/Nrep * log[P(spike|Stim)] + NsilencesperStim/Nrep * log[P(no spike|Stim)])

% Author: Jeff Fitzgerald, August 2010 (updated on Sept 5, 2010)


function f = logloss(p, stim, resp, order)

[Nsamples,Ndim] = size(stim);

ptemp = p(2:Ndim+1);
if order>1
    J = reshape(p(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim])';
end

if order==1
    f1 = 1+exp(p(1)+stim*ptemp');
    f0 = 1+exp(-p(1)-stim*ptemp');
else
    f1 = 1+exp(p(1)+stim*ptemp'+sum(stim.*(stim*J),2));
    f0 = 1+exp(-p(1)-stim*ptemp'-sum(stim.*(stim*J),2));
end

NspikesperStim = resp;  % Nsamples x 1
F1 = NspikesperStim.*log(f1);
F0 = (1-NspikesperStim).*log(f0);
F1(isnan(F1)) = 0;
F0(isnan(F0)) = 0;
f = mean(F0+F1);