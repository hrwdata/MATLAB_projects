% % %                                                   % % % 
% % %   PRIMAL - DUAL SIMPLEX METHOD                    % % % 
% % %                                                   % % % 
% % %   min z = cTx, s.t. Ax = b, x >= 0                % % % 
% % %       with an initial dual feasible               % % %
% % %       solution available                          % % %  
% % %                                                   % % %  
% % %   lam - dual feasible solution                    % % % 
function [x_star z_star lam_star] = primal_dual(A,b,c,lam)
    
% step 0

    data = [A b; c' 0];
    
    [m n] = size(A);
    
    std_ratio = 0;
    
% step 1
    
    % feasible lam is available
    
% step 2 

    % set P = P_l = { i: c_l = c'-lam'A = 0}

    c_l = c'-lam'*A;
    
    c_l = [c_l nan(1,m)];
    
    P_l = ~(c_l~=0);
    
    % print P, P_c
    P = P_l(1:n);
    
    P_c = ~P;
    
    P_c = ~P_l;
    
    
    % get ARP
        
    A_ARP = [A eye(m)];
        
    B_ARP = eye(m); % initial B^-1 = I
    
    %initial c_ARP
    c_ARP = [zeros(n,1); ones(m,1)]; 
    
    %basic entries of c
    c_B = ones(m,1);
    
    % u=c_B'B^-1 cf. Section 4.3
    u = c_B'*inv(B_ARP(:,1:m));
    
    % initial r_ARP
    r_ARP = zeros(m+n,1);
    
    D_ARP = A_ARP(:,P_c);
    
    r_ARP(P_c) = (c_ARP(P_c)'-u*D_ARP)';
                    
    % ARP objective function
    z_ARP = c_B'*inv(B_ARP(:,1:m))*b;

    % print ARP tableau
    ARP = [A_ARP b;r_ARP' -z_ARP];
        
    iter=1;
    
    iter_3=1;
    
    
% step 3

    r_ARP_0 = r_ARP;
    
    r_ARP_0(n+1:end)=NaN(m,1);
    
    for(iter=1:1000)

    % solve ARP
    r_ARP_P = r_ARP(1:n);
    r_ARP_P(~P) = NaN;
    r_ARP_P;
    
    condition_3 = all(r_ARP_P(P) >= 0);
    
    P_1 = P;
    
    if ~condition_3
        [~, q] = min(r_ARP_P);
    end
    
    if isnan(std_ratio)
        warning('primal_dual:InitialDualInfeasible', 'Initial dual solution is infeasible.');
        break
    end  
    
    while(condition_3~=1&&iter_3<=1000)
        
        % standard simplex method ratio
        std_ratio = ARP(1:m,end)./ARP(1:m,q);
        
        std_ratio(ARP(1:m,q)<0) = NaN;
        
        [std_ratio p] = min(std_ratio);
        
        if isnan(std_ratio)
            warning('primal_dual:Unbounded', 'Unbounded.');
            break
        end    
        
        % pivot about (p,q) 
        % row operations
        ARP(p,:) = ARP(p,:)/ARP(p,q);
        
        for index=[1:p-1 p+1:m+1]
        
            ARP(index,:) = ARP(index,:)- ARP(p,:).*ARP(index,q);
        
        end
            
        % new r_ARP
        r_ARP = ARP(end,1:m+n)';
            
        % step 3 simplex iteration
        iter_3=iter_3+1;
            
        P_1(p)=0;
        
        r_ARP_3 = r_ARP(1:n);
        
        r_ARP_3(~P_1)=NaN;
        
        [~, q] = min(r_ARP_3);
        
        % all r>0?
        condition_3 = all(r_ARP_3(P_1)>=0);
            
        % print ARP
        ARP;
            
        if(iter_3==10000)
            
            warning('primal_dual:Step3Timeout', 'Time expired: 10,000 iterations at Step 3.');
            
            break
        
        end
        
    end

% step 4
    
    % basic variable index
    bvi = zeros(1,m+n);
    
    for jndex = 1:m+n
        
        [bvi_1,~,bvi_0] = find(ARP(:,jndex));
    
        bvi_0=[bvi_0; size(bvi_1,1)];
        
        bvi(jndex) = sum(all(bvi_0==1)*bvi_1);
    
    end
    
    % x,y
    x = zeros(1,m+n);
    
    bvi=bvi;
    
    ARP_1 = ARP(1:m,end)';
    
    x(logical(bvi)) =  ARP_1(bvi(logical(bvi)));
    
    bvi_2 = bvi;
    
    bvi_2(bvi_2==0)=[];
    
    y = x(n+1:end);
    
    y_0 = all(y==0);
    
    condition_4 = y_0;
    % y=0? If no, continue.
    
    % if all y=0 -> optimal solution = x
    if(y_0~=0) 
            
        x_star=x(1:n)';

        z_star=c'*x_star;

        lam_star=lam;

        break 
     
    end
         
    % step 5    
               
    r_ARP_0 = all(r_ARP>=0);
    
    r_ARP;
    
    condition_5 = r_ARP_0;
            
    % r_ARP>=0?     % If no, continue.
    if(r_ARP_0~=0) 
        
        warning('primal_dual:PrimalInfeasible', '(P) is infeasible.');
        
        break 
    
    end
            
    
% step 6

    P=P;
    
    c_ARP_1=c_ARP;
    
    c_ARP_1(bvi(logical(bvi))) =  c_ARP(logical(bvi));
    
    
    % [u = (cB_ARP)B^-1 cf. Section 4.3]
    c_ARP_2=c_ARP_1(bvi(logical(bvi)))';
    ARP_2 = ARP(1:m,n+1:end-1);
    u = c_ARP_2(1:m)*ARP_2;
    
    % set up for c'-lam'*A ratio
    r_ARP_r = r_ARP;
    
    r_ARP_r(n+1:end)=NaN(m,1);
    
    % (c'-lam'*A)/(-r_ARP)'
    
    r_ARP_r(r_ARP_r>=0)=NaN;
    
    ratio_l = c_l./-r_ARP_r';
    
    ratio_l(P_l) = NaN; 
    
    % min { (cj - lam'*aj)/(- rj_ARP): rj_ARP < 0, j in P_c}
    [eps q] = min(ratio_l(1:n));
    
    % eps = 0 if eps = emptyset
    if isnan(eps)
        warning('primal_dual:Infeasible', 'Infeasible.');
        break        
    end
        
    % new lam = old lam + eps*u
    lam = lam + eps*u';
    
    % step 2 
        % set P = P_l = { i: c_l = c'-lam'A = 0}
    r_ARP_2 = r_ARP;
    
    r_ARP_2(n+1:end)=NaN;
        
    c_l = c_l+eps*r_ARP_2';
        
    P_l = ~(c_l~=0);
    
    P_l(n+1:end) = 0;
    
    % print P,P_c
    P = P_l(1:n);
    
    P_c = ~P;
    
    P_c = ~P_l;
    
    % get ARP
        
    A_ARP = ARP(1:m,:);
       
    B_ARP = A_ARP(1:m,n+1:end-1);
                
    c_B = c_ARP(logical(bvi));
        
    D_ARP = A_ARP(:,~logical(bvi));
    
    z_ARP = -ARP(m+1,m+n+1);
    
    b_ARP = ARP(:,m+n+1);
    
    iter=iter+1;
    
    ARP;
    
    if(iter==10000)
        
        warning('primal_dual:Step6Timeout', 'Time expired: 10,000 iterations at Step 6.');
        
        break
        
    end
    
    end
