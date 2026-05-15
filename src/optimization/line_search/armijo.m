%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%        Armijo's Line Search Algorithm     %%%
%%%          (P) minimize f(x+alpha*d)        %%%
%%%    arguments                              %%%
%%%    f  - symbolic scalar function          %%%
%%%    x0 - numeric initial point             %%%
%%%    d  - search direction                  %%%
%%%    epsilon - bound for gradient           %%%
%%%    eta - factor to increase alpha         %%%
%%%                                           %%%
%%%    Author: Harrison Watts                 %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = armijo(f,x0,d,epsilon,eta)
    
    % epsilon = 1/2
    if nargin < 4 || isempty(epsilon)
        epsilon = 1/2;
    end
    
    % eta = 2
    if nargin < 5 || isempty(eta)
        eta = 2;
    end
    
    x0 = x0(:);
    d = d(:);
    vars = symvar(f);
    grad = gradient(f, vars);
    
    % initial alpha
    alpha = 1;

    % starting objective value
    phi_0 = evaluate_symbolic_at_point(f, vars, x0);
    
    % alpha derivative
    g0 = evaluate_symbolic_at_point(grad, vars, x0);
    phi_prime_0 = g0(:)'*d;
    
    % while condition is not met
    condition_met = 0;
    
    while ~condition_met
        
        % new location
        x_next = x0 + alpha*d;
        
        % objective function value
        phi = evaluate_symbolic_at_point(f, vars, x_next);
        
        % Armijo's rule
        armijos_rule = double(phi_0 + epsilon*phi_prime_0*alpha);

        condition_met = (phi <= armijos_rule);
        
        % alpha is too large
        if ~condition_met
            
            alpha = alpha/3;

            % objective function value at updated alpha
            x_next = x0 + alpha*d;
            phi = evaluate_symbolic_at_point(f, vars, x_next);
            
            % Armijo's rule
            armijos_rule = double(phi_0 + epsilon*phi_prime_0*alpha);
            
            condition_met = (phi <= armijos_rule);
            
        % alpha is too small
        else 
            
            alpha = alpha*2;

            % objective function value at updated alpha
            x_next = x0 + alpha*d;
            
            % curvature condition
            phi_eta = evaluate_symbolic_at_point(f, vars, eta*x_next);
            
            curvature_condition = double(phi_0+epsilon*phi_prime_0*eta*alpha);
            
            condition_met = (phi_eta > curvature_condition);

        end
        
    end
    
    alpha = alpha;
   
end
