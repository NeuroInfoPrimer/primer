 
function [df] = df1dim(x, stim, avgs, order)    
% Must accompany linmin. 
 
global pcom; % Defned in linmin. 
global xicom; 
global nrdfun; 

df  = dot(feval(nrdfun, pcom + x.*xicom, stim, avgs, order),xicom); 