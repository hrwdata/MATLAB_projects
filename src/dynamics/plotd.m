function X = plotd(yprime,a,b,c,d,h)
% Direction fields and solutions using quiver
% yprime = dy/dx 
% [a,b] = x interval
% [c,d] = y interval
% h = stepsize
    % optional
    if nargin < 6 || isempty(h)
        h = (b-a)./10;
    end

    % create a grid of points where slope vectors will be drawn
    [t,y]=meshgrid(a:h:b,c:h:d);
    
    % dimensions of grid 
    m=size(t);
    n=size(y);
    
    % differential equation
    dt = ones(m);
    dy = yprime(t,y);
    
    % make arrows same size
    veclength = sqrt(dt.^2+dy.^2);
    dt = dt./veclength;
    dy = dy./veclength;
    
    % plot slope field
    X = quiver(t,y,dt,dy);
end
