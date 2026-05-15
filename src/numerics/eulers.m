function [x y f table] = eulers(yprime,a,b,y0,h)
% Euler's method on [a,b]
% Author: Harrison Watts
% Date published: Dec. 18, 2021
% 
% y' = f(x,y), y(a)=y0
% yprime = @(x,y) function that describes y'
% a = x0, starting point for x
% b = endpoint of interval of evaluation
% y0 = starting point for y
% h = stepsize
% 
% x = x0 + nh
% y = solution values
% f = slopes at each iteration
% table = printed table of results
% 

% initial values
x = a:h:b;
y = y0;
f = yprime(x(1),y(1));

% number of steps
N = size(x,2); %size(x,2) gives #columns in x
% start loop
for n=1:N-1

    % Euler's method 
    % y(n+1) = y(n) + h*f(x(n),y(n))
    y(n+1) = y(n)+h*f(n);

    % f(n+1) = f(x(n+1),y(n+1))
    f(n+1) = yprime(x(n+1),y(n+1));
end
    
% print results
table_1 = [x; y; f]';
rowNames = arrayfun(@num2str,0:N-1,'Uni',0);
colNames = {'x','y','f'};
table = array2table(table_1,'RowNames',rowNames,'VariableNames',colNames);    
end
