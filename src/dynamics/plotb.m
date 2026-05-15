function plotb(pprime,p0,h1,h2,hinc)
%% Bifurcation Diagram
% pprime, p0 = differential equation
% h1 = start, h2 = end, hinc = increment

src_dir = fileparts(fileparts(mfilename('fullpath')));
numerics_dir = fullfile(src_dir, 'numerics');
original_dir = pwd;
cleanup_obj = onCleanup(@() cd(original_dir));
cd(numerics_dir);

figure;
hold on

for h = h1:hinc:h2
    
    % perform 230 iterations of Euler's method
    b = 230*h;
    [~, p] = eulers(pprime,0,b,p0,h);
    
    % plot 201-230
    h_time = repmat(h,30,1);
    plot(h_time,p(202:231),'.','Color','blue')
end

xlabel('h')
ylabel('limit values')
xlim([h1 h2])
hold off
end