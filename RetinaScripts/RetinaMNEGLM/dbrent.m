%| FUNCTION: dbrent 
%| 
%| PURPOSE: Given a function f and its derivative function df, and  
%|   given a bracketing triplet of abscissas ax, bx, cx [such that  
%|   bx is between ax and cx, and f(bx) is less than both f(ax) and  
%|   f(cx)], this routine isolates the minimum to a fractional precision  
%|   of about tol using a modification of Brent's method that uses  
%|   derivatives. The abscissa of the minimum is returned as xmin, and 
%|   the minimum function value is returned as rval, the returned  
%|   function value. 
%| 
%| REFERENCE:  Numerical recipes in C 
%| 
 
 
function [rval, xmin] = dbrent(ax, bx, cx, f, df, tol, stim, resp, order, avgs) 
 
ITMAX = 50; 
ZEPS  = 1.0e-10; 
e     = 0.0; 
 
if ax < cx, 
   a = ax; 
   b = cx; 
else  
   a = cx; 
   b = ax; 
end 
 
x  = bx;          w  = bx; v  = bx; 
fx = feval(f,x, stim, resp, order);  fw = fx; fv = fx; 
dx = feval(df,x, stim, avgs, order); dw = dx; dv = dx; 
%xi=dx;
 
for iter=1:ITMAX, 
   xm   = 0.5*(a+b); 
   tol1 = tol*abs(x)+ZEPS; 
   tol2 = 2.0*tol1; 
   if abs(x-xm) <= (tol2-0.5*(b-a)) 
      xmin  = x; 
      rval  = fx; 
      return; 
   end 
 
   if (abs(e) > tol1)  
      d1 = 2.0*(b-a); % Initialize these d's to an out-of-bracket value 
      d2 = d1; 
      if (dw ~= dx), 
         d1 = (w-x)*dx/(dx-dw); % Secant method with one point. 
      end 
      if (dv ~= dx), 
         d2 = (v-x)*dx/(dx-dv); % Secant method with the other point. 
      end 
       
      %|  Which of these two estimates of d shall we take? We will  
      %|  insist that they be within the bracket, and on the side  
      %|  pointed to by the derivative at x: 
      u1   = x+d1; 
      u2   = x+d2; 
      ok1  = ((a-u1)*(u1-b) > 0.0) & (dx*d1 <= 0.0); 
      ok2  = ((a-u2)*(u2-b) > 0.0) & (dx*d2 <= 0.0); 
      olde = e; % Movement on the step before last. 
      e    = d; 
       
      %|  Take only an acceptable d, and if both are acceptable,  
      %|  then take the smallest one. 
      if ok1 | ok2, 
         if ok1 & ok2, 
            if abs(d1) < abs(d2), d = d1; else d = d2; end   
         elseif ok1, 
            d = d1; 
         else 
            d = d2; 
         end; 
          
         if abs(d) <= abs(0.5*olde),  
            u=x+d; 
            if (u-a < tol2) | (b-u < tol2), 
               d = tol1*signNR(xm-x); 
            end 
         else   % Bisect, not golden section. 
            %|  Decide which segment by the signNR of the derivative. 
            if dx >= 0, e=a-x; else e=b-x; end 
            d = 0.5*e; 
         end 
      else 
         if dx >= 0, e=a-x; else e=b-x; end 
         d = 0.5*e; 
      end 
   else 
      if dx >= 0, e=a-x; else e=b-x; end 
      d = 0.5*e; 
   end 
   if abs(d) >= tol1, 
      u  = x + d; 
      fu = feval(f,u, stim, resp, order); 
   else 
      u  = x + tol1*signNR(d); 
      fu = feval(f,u, stim, resp, order); 
      %|  If the minimum step in the downhill direction takes us  
      %|  uphill, then we are done. 
      if fu > fx,   
         xmin = x; 
         rval = fx; 
         return; 
      end 
   end 
    
   %|  Housekeeping 
   du = feval(df, u, stim, avgs, order); 
   if fu <= fx, 
      if u >= x, a=x; else b=x; end; 
      v = w; fv = fw; dv = dw; 
      w = x; fw = fx; dw = dx; 
      x = u; fx = fu; dx = du; 
   else 
      if u < x, a=u; else b=u; end; 
      if (fu <= fw) | (w == x), 
         v = w; fv = fw; dv = dw; 
         w = u; fw = fu; dw = du; 
      elseif (fu < fv | v == x | v == w), 
         v = u; fv = fu; dv = du; 
      end; 
   end; 
end; 
 
disp('Too many iterations in routine dbrent'); 
 
xmin = x; 
rval = fx;  % hopefully never get here.