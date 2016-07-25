function [ax, bx, cx, fa, fb, fc]=mnbrak(ax,bx,f, stim, resp, order, avgs)
%brackets a minimum of f
%must be supplied with f=@(x) 

% Constants
GOLD=1.618034;
GLIMIT=100.0;
TINY=1.0e-20;

fa=f(ax, stim, resp, order);
fb=f(bx, stim, resp, order);

if fb>fa,
    dum=ax ; ax=bx ; bx=dum ;
    dum=fb ; fb=fa ; fa=dum ;
end;

cx=bx+GOLD*(bx-ax);
fc=f(cx, stim, resp, order);

while fb>fc
    r=(bx-ax)*(fb-fc);
    q=(bx-cx)*(fb-fa);
    u= bx -((bx-cx)*q-(bx-ax)*r)/(2.0*max([abs(q-r),TINY])*sign(q-r));
    ulim=(bx)+GLIMIT*(cx-bx);
    if (bx-u)*(u-cx)>0.0,
        fu=f(u, stim, resp, order) ;
        if (fu < fc)
            ax=bx;
            bx=u;
            fa=fb;
            fb=fu;
            return;
        elseif (fu> fb)
            cx=u;
            fc=fu;
            return;
        end
        u=cx+GOLD*(cx-bx);
        fu=f(u, stim, resp, order);
    elseif (cx-u)*(u-ulim)>0.0,
        fu=f(u, stim, resp, order) ;
        if fu<fc
            bx=cx;
            cx=u;
            u=u+GOLD*(u-cx);
            fb=fc;
            fc=fu;
            fu=f(u, stim, resp, order);
        end
    elseif ((u-ulim)*(ulim-cx))>=0.0
        u=ulim;
        fu=f(u, stim, resp, order);
    else
        u=(cx)+GOLD*(cx-bx);
        fu=f(u, stim, resp, order) ;
    end ;
    
    ax=bx ; bx=cx ; cx=u ;
    fa=fb ; fb=fc ; fc=fu ;
end

