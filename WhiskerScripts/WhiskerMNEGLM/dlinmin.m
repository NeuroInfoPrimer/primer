%| FUNCTION: dlinmin 
%| 
%| PURPOSE:  Given an n-dimensional point p[1..n] and an  
%|   n-dimensional direction xi[1..n], moves and resets p to where  
%|   the function func(p) takes on a minimum along the direction xi  
%|   from p, and replaces xi by the actual vector displacement that  
%|   p was moved. Also returns as fret the value of func at the returned  
%|   location p. This is actually all accomplished by calling the 
%|   routines mnbrak and brent. 
%| 
%| REFERENCE:  Numerical recipes in C 
%| 
 
function [p, xi, fret] = dlinmin(p, xi, func, dfunc, stim, resp, order, avgs)
 
TOL    = 2.0e-4;  % Tolerance passed to brent. 

global pcom xicom nrfunc nrdfun; 
nrfunc = func; 
nrdfun = dfunc; 
pcom   = p; 
xicom  = xi; 
 
ax     = 0.0;  % Initial guess for brackets. 
xx     = .2;%2*rand();%2.0; 

ftemp=@f1dim;
ftemp2=@df1dim;

[ax, xx, bx, fa1, fc1, fb1] = mnbrak(ax, xx, ftemp, stim, resp, order, avgs);
%[ax,xx,bx]
%plotalonglinemin(ax,bx,ftemp);
[fret, xmin] = dbrent(ax,xx,bx,ftemp,ftemp2,TOL, stim, resp, order, avgs); 
%xmin
%plot([xmin],[fret],'ro','LineWidth',7); 
%drawnow;
%hold off;
xi     = xi.*xmin;
p      = p + xi;

