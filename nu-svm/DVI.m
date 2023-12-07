function [R1,R2,delta1,t1,t2,Emax,Ebar_i0,Ebar_j0,D,length_delta1_Nbar]= DVI(kernel_type,s,delta0,Q,alpha0,nu1,traY)
l=size(Q,1); 
C1=1/l;
N_bar = find(delta0< -alpha0 | delta0 > C1-alpha0); 
delta1 = delta0; 
lambda = 0.5; 
solve_delta_times = 4; 
length_delta1_Nbar=0;

tic
if strcmp(kernel_type,'rbf')
    Qalpha = Q*alpha0;
    P = 2*Q;  
    p = 4*Qalpha;  
else
    [Qii,i0] = max(diag(Q));
    P = 0.5*sqrt(Qii)*Q;  
    p = sqrt(Qii)*Q* alpha0+1/2*Q(:,i0);  
end 

if isempty(N_bar)
    solve_time = 0;
else
    if s~=2 && (length(N_bar)< lambda*l )  
        delta1_Nbar = []; 
        N = 1:l;
        N = N(ismember(N,N_bar)==0); 
        solve_time = 0;
        while length(N_bar) < lambda*l || isempty(delta1_Nbar)==1 
            n = min(solve_delta_times*length(N_bar),lambda*l);
            sample = randperm(length(N),round(n));
            index = N(sample);
            N_bar = [N_bar;index'];  
            N(sample) = [];  
            delta0_N = delta0;
            delta0_N(N_bar) = []; 
            
            P_Nbar_Nar = P(N_bar,N_bar); 
            P_Nbar_N = P(N_bar,N);
            p_bar = P_Nbar_N*delta0_N + p(N_bar);
            A = -ones(1,length(N_bar));
            b = -nu1+sum(delta0_N)+sum(alpha0);
            Aeq = [];
            beq = [];
            lb = -alpha0(N_bar);
            ub = C1-alpha0(N_bar);
            x0=[];
            delta1_Nbar = quadprog(P_Nbar_Nar,p_bar,A,b,Aeq,beq,lb,ub,x0);
            length_delta1_Nbar=length(delta1_Nbar);
            solve_time = solve_time +1;
            if isempty(delta1_Nbar)==0 
                delta1(N_bar)=delta1_Nbar; 
                break;
            end
            if solve_time>=6 
                break;
            end
        end
    else   
        A = -ones(1,l);
        b = -nu1+sum(alpha0);
        Aeq = [];
        beq = [];
        lb = -alpha0;
        ub = C1-alpha0;
        [delta1]=quadprog(P,p,A,b,Aeq,beq,lb,ub);
        length_delta1_Nbar=l;
        solve_time = 1;
    end
end
t1= toc;

tic

i0 = floor(l-nu1*l);
j0 = ceil(l-nu1*l);
if strcmp(kernel_type,'rbf')
    E = Q*(alpha0+0.5*delta1);
    [Ebar,I] = sort( E , 'descend' );
    Emax = Ebar(1);
    Ebar_i0 = Ebar(i0);
    Ebar_j0 = Ebar(j0);
    D = 2*sqrt( abs( (alpha0+0.25*delta1)'*Q*delta1 ) );
    R1 =  E > ( E(I(i0))+D);
    R2 =  E < ( E(I(j0))-D);
else
    s
    size(alpha0)
    size(Q)
    size(delta1) 
    
    theta=(Q*alpha0)'*delta1+1/4*delta1'*Q*delta1; 
    norm_phix=sqrt(diag(Q));
    y_w_phi_down=Q*(alpha0+0.5*delta1)-sqrt(abs(theta))*norm_phix;
    y_w_phi_up=Q*(alpha0+0.5*delta1)+sqrt(abs(theta))*norm_phix;
    [y_w_phi_down_descend,~]=sort(y_w_phi_down,'descend');
    [y_w_phi_up_descend,~]=sort(y_w_phi_up,'descend');
    d_up_down=y_w_phi_up_descend(i0);
    d_down_up=y_w_phi_down_descend(j0);
    R1=y_w_phi_down>d_up_down;
    R2=y_w_phi_up<d_down_up;
    
    Emax = y_w_phi_down_descend(1);
    Ebar_i0 = d_up_down;
    Ebar_j0 = d_down_up;
    D = theta;
end
t2 = toc;

