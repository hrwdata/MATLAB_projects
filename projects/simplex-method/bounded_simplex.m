% % %                                                       % % % 
% % %   BOUNDED-SOLUTION SIMPLEX METHOD                     % % % 
% % %                                                       % % % 
% % %   min z = cTx, s.t. Ax = b, d<= x <= h                % % % 
% % %   with and without an initial                         % % %  
% % %   solution available                                  % % %  
% % %                                                       % % % 
% % %                                                       % % % 
% % %   lam - dual feasible solution                        % % % 
function [x_star z_star] = bounded_simplex(A,b,c,d,h)

 b = b-A*d
 
 h = h-d
    % Phase I
    
    [m n] = size(A);
    
    c_I = [zeros(1,n) ones(1,m) 0];
    % c_I = [c' zeros(1,m) 7]
    
    tableau = [A eye(m) b; c_I];
    
    tableau(end,:) = tableau(end,:)-sum(tableau(1:m,:),1)
    
    
    % basic variable index
    bvi = zeros(1,m+n);
    
    for jndex = 1:m+n
        
        [bvi_1,~,bvi_0] = find(tableau(:,jndex));
    
        bvi_0=[bvi_0; size(bvi_1,1)];
        
        bvi(jndex) = sum(all(bvi_0==1)*bvi_1);
    
    end
    
    q_1=[];
            X_1 = zeros(1,m+n);

    for iter=1:1000
    

    X_2 =  tableau(1:m,end)';
    
    bvi_1 = bvi;
    
    bvi_1(bvi==0) = [];
    
    X_1(logical(bvi)) =  X_2(bvi_1);
    
    % pivot column q
    
    r_1 = tableau(end,1:m+n);
    
    r_1(logical((tableau(end,1:m+n)>=0)))=NaN;
    
    if ~isempty(q_1)
        r_1(q_1)=NaN;
    end
    
    [~, q] = min(r_1)
    
    r_2=r_1;
    
    r_2(isnan(r_1))=[];
    
    if isempty(r_2)
        break
    end
    
    % pivot row p
    
    eps = NaN(3,2);
    
    h_I = [h; Inf(m,1)];
    
    eps(1,1) = h(q);
    
    ratio_1 = tableau(1:m,end)./tableau(1:m,q);
    
    if ~isempty(min(ratio_1(tableau(1:m,q)>=0)))
        
        ratio_1(tableau(1:m,q)<0)=NaN;

        [eps(3,1) eps(3,2)] = min(ratio_1);
    
    end
    
    ratio_2=(h(q)-b)./tableau(1:m,q);
    
    if ~isempty(min(ratio_2(tableau(1:m,q)>0)))
        
        ratio_2(tableau(1:m,q)<=0)=NaN;
    
        [eps(2,1) eps(2,2)] = min(ratio_2);
    
    end
    
    eps(eps<=0)=NaN;
            
    [eps_4 p_4] = min(eps(:,1))
    
    p = eps(p_4,2)

    
    % bound reached
    if p_4==1
        
        X_1(q)=eps_4;
        
        % fix nonbasic
        q_1=[q_1 q];
        
        tableau
        
    end
    
    
    % adjust bounds then pivot
    if p_4==2
        
        X_1(q)=eps_4;
        
        tableau(:,end)=tableau(:,end)-X_1(q)*tableau(:,q);
        
        tableau(:,q)=tableau(:,q)-X_1(q)*tableau(:,q);
       
            
        tableau
        
        q_1 =[q_1 q];
        
    end
    
    % standard simplex, bound not reached
    if p_4==3
    
        % pivot about (p,q) 
        % row operations
        tableau(p,:) = tableau(p,:)/tableau(p,q);
        
        for index=[1:p-1 p+1:m+1]
        
            tableau(index,:) = tableau(index,:)- tableau(p,:).*tableau(index,q);
        
        end
            
        tableau
        
    end
    
    if all(tableau(end,1:end-1) >= 0)
        if tableau(end,end)~=0
            'infeasible'
            break
        end
        break
    end
    
    end
    
    % Phase II
    tableau(:,n+1:end-1)=[]
        
    % basic variable index
    bvi = zeros(1,n);
    
    for jndex = 1:n
        
        [bvi_1,~,bvi_0] = find(tableau(:,jndex));
    
        bvi_0=[bvi_0; size(bvi_1,1)];
        
        bvi(jndex) = sum(all(bvi_0==1)*bvi_1);
    
    end
    
    bvi
    
    c_II = [c' 0];
    % c_I = [c' zeros(1,m) 0]
    
    
    c_II_1 = c';
    
    c_II_1(~logical(bvi))=0;
    
    bvi_1=bvi(bvi~=0);
    
    c_II_1(bvi_1)=c_II_1(logical(bvi));
    
    tableau(end,:)=c_II-c_II_1(1:m)*tableau(1:m,:)

    
    % basic variable index     
    bvi = zeros(1,n);
    
    X_1 = zeros(1,n);
    
    for iter=1:5
        
        iter
        
        for jndex = 1:n
        
            [bvi_1,~,bvi_0] = find(tableau(:,jndex));
    
            bvi_0=[bvi_0; size(bvi_1,1)];
        
            bvi(jndex) = sum(all(bvi_0==1)*bvi_1);
    
        end
         
        if all(tableau(end,1:end-1) >= 0)
            x_star = X_1'+d
            z_star = c'*x_star
            break
        end
        
    
        X_2 =  tableau(1:m,end)';
    
        bvi_1 = bvi;
    
        bvi_1(bvi==0) = [];
    
        X_1(logical(bvi)) =  X_2(bvi_1);
    
        X_1(q_1)=h(q_1);
    
    
    
        % pivot column q
    
        r_1 = tableau(end,1:n);
    
        r_1(logical((tableau(end,1:n)>=0)))=NaN;
    
        [~, q] = min(r_1)
    
        r_2=r_1;
    
        r_2(isnan(r_1))=[];
    
        if isempty(r_2)
            x_star = X_1'+d
            z_star = c'*x_star
            break
        end
    
        % pivot row p
    
        eps = NaN(3,2);
    
        h_I = [h; Inf(m,1)];
    
        eps(1,1) = h(q);
    
        ratio_1 = tableau(1:m,end)./tableau(1:m,q);
    
        if ~isempty(min(ratio_1(tableau(1:m,q)>=0)))
    
            [eps(3,1) eps(3,2)] = min(ratio_1(tableau(1:m,q)>=0));
    
        end
     
        ratio_2=(h(q)-b)./tableau(1:m,q);
    
    
        if ~isempty(min(ratio_2(tableau(1:m,q)>0)))
            ratio_2(tableau(1:m,q)<=0)=NaN;
            if ~isempty(ratio_2(ratio_2<=0))
                ratio_2(ratio_2<=0)=NaN;
                [eps(2,1) eps(2,2)] = min(ratio_2);
            end
        end
            
        [eps_4 p_4] = min(eps(:,1));
    
        p = eps(p_4,2)
    
        % bound reached
        if p_4==1
        
            X_1(q)=eps_4;
        
            tableau(:,end)=tableau(:,end)-X_1(q)*tableau(:,q);
        
            tableau(:,q)=tableau(:,q)-X_1(q)*tableau(:,q);
        
            % fix nonbasic
            q_1=q;
        
            tableau
        
        end
    
    
        % adjust bounds then pivot
        if p_4==2
        
            X_1(q)=eps_4;
        
            tableau(:,end)=tableau(:,end)-X_1(q)*tableau(:,q);
        
            tableau(:,q)=tableau(:,q)-X_1(q)*tableau(:,q);
        
            % pivot about (p,q) 
            % row operations
            tableau(p,:) = tableau(p,:)/tableau(p,q);
        
            for index=[1:p-1 p+1:m+1]
        
                tableau(index,:) = tableau(index,:)- tableau(p,:).*tableau(index,q);
        
            end
            
            tableau
        
        end
    
        % standard simplex, bound not reached
        if p_4==3
    
            % pivot about (p,q) 
            % row operations
            tableau(p,:) = tableau(p,:)/tableau(p,q);
        
            for index=[1:p-1 p+1:m+1]

            tableau(index,:) = tableau(index,:)- tableau(p,:).*tableau(index,q);

            end
            tableau
        
        end
    
    end
   
end