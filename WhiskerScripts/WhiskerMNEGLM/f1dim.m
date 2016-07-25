% Move pcom (x units in xicom direction), and then evaluate the function there. 
 
function [f] = f1dim(x, stim, resp, order)    
% Must accompany linmin. 

global pcom; % Defned in linmin. 
global xicom; 
global nrfunc; 

f  = feval(nrfunc, pcom + x.*xicom, stim, resp, order); 
