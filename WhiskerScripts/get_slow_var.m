%
% Use the phase (p) to find the turning points of the whisks (tops and
% bottoms), and calculate a value on each consecutive whisk using the 
% function handle (operation).  The values are calculated twice per 
% whisk cycle using both bottom-to-bottom and top-to-top.  The values
% are linearly interpolated between whisks.
%
function [out, tops, bottoms] = get_slow_var( sig, p, operation )

    % find crossings
    tops = find(p(1:end-1)<0 & p(2:end)>=0);
    bottoms = find(p(1:end-1)>=pi/2 & p(2:end)<=-pi/2);

    out = zeros( [1 length(sig)] );

    % evaluate at transitions
    temp = [];
    pos = [];
    for j = 2:length(tops)
        vals =  sig( tops(j-1):tops(j) );
        temp(end+1)  =  operation(vals);
    end
    if length(tops) > 1
     pos =    round( tops(1:end-1) + diff(tops)/2);
    end
    for j = 2:length(bottoms)
        vals =  sig( bottoms(j-1):bottoms(j) );
        temp(end+1)  =  operation(vals);
    end
    if length(bottoms) > 1
     pos   = [pos round(bottoms(1:end-1) + diff(bottoms)/2)];
    end

    % sort everything
     [pos,i] = sort(pos);
     pos = [1 pos length(sig)];

     if isempty(temp)
         temp = operation( sig ) * [1 1] ;
     else
       temp = [temp(i(1)) temp(i) temp(i(end))];
     end

    % make piecewise linear signal
    for j = 2:length(pos)
        in = pos(j-1):pos(j);
        out(in) = linspace( temp(j-1), temp(j), length(in) );
    end

