clear;
% close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Homothetic Tube-based EMPC algorithm %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    format long

    global t N V A B zs Omegavertex coeff;
    
    % system description
    A = [0.7776 -0.0045; 26.6185 1.8555];
    B = [-0.0004; 0.2907];
    Bw = [-0.0002 0.0893; 0.139 1.2267];
    C = [0 1];
    
    load('Sets.mat')
    load('StorageCoeff.mat');
    
    nx = 2;
    nu = 1;
    ny = 1;
    
    % constraints
    xh = [0.5; 5];
    xl = [-0.5; -5];
    
    uh = 15;
    ul = -15;
    
    wh = [0.2; 0.05];
    wl = [-0.2; -0.05];
    
    etah = 0.01;
    etal = -0.01;
    
    Aw = [eye(2); -eye(2)];
    bw = [wh; -wl];
    Pw = Polyhedron(Aw,bw);
    
    % steady state   
    load('SteadyState.mat')
    
    % tightened constranits
    zh = [0.4811; 4.139];
    zl = [-0.4811; -4.139];
    
    vh = 9.997;
    vl = -9.997;
    
    % horizons
    M = 100;
    N = 20;
    
    
    % state disturbance
    load('Noises.mat');
%     W1 = wl(1) + (wh(1)-wl(1))*rand(1,M);
%     W2 = wl(2) + (wh(2)-wl(2))*rand(1,M);
%     W = [W1; W2];
%     Eta = etal + (etah-etal)*rand(1,M+1);

%     W = kron(ones(1,M),wh);
%     Eta = etal*ones(1,M+1);
    
    
%     W = zeros(2,M);
%     Eta = zeros(1,M+1);
    
    % input and output data
    V = zeros(nu,M);
    U = zeros(nu,M);
    Gamma = zeros(1,M);
    Y = zeros(ny,M+1);
    VavSum = zeros(M,1);
    AsyAve = zeros(M,1);
    VavSumUpper = zeros(M,1);
    AsyAveUpper = zeros(M,1);
%     VavSumLower = zeros(M,1);
%     AsyAveLower = zeros(M,1);
    
    % initial state
    X = [[-0.2; 2.5] zeros(nx,M)];
    Xesti = [X(:,1)+[0.005; 0.1] zeros(nx,M)];
    Z = zeros(nx,M);
    Y(1) = C*X(:,1) + Eta(1);
    
    
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10);
%     options = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10);

    % average time for +,*
    niter = 2000;
    a = rand(100);
    b = rand(100);
    tt = [];
    tic; for i=1:niter a+4; end; tt(1)=toc/niter;
    tic; for i=1:niter a+b; end; tt(2)=toc/niter;
    tic; for i=1:niter a.*b; end; tt(3)=toc/niter;
    t_arith = mean(tt);
    t_flop = t_arith/prod(size(a));
        
    t_fmincon = zeros(M,1);
    flops_fmincon = zeros(M,1);


    t = 0;
    
    % bound
    lb = [kron(ones(t+N+1,1),zl); kron(ones(N,1),vl); kron(ones(N+1,1),0)];
    ub = [kron(ones(t+N+1,1),zh); kron(ones(N,1),vh); kron(ones(N+1,1),1)];
    
    % initial state constraint
    AEMPCineq = [-Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N+N+1); Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N+N+1)];
    bEMPCineq = [Xi.H(:,nx+1)-Xi.H(:,1:nx)*Xesti(:,t+1); Xi.H(:,nx+1)+Xi.H(:,1:nx)*Xesti(:,t+1)];
    
    
    % terminal equality constraint
    AEMPCeq = [zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N+N+1)];
    bEMPCeq = zs;
    
    % reachable set
    Omegavertex = Omega.V';
    Rset = Xesti(:,t+1) + E;
    
    Rset_A = [Rset.H(:,1:nx); C; -C];
    Rset_b = [Rset.H(:,nx+1); Y(t+1)+etah; -(Y(t+1)+etal)];
        
    Rset = Polyhedron(Rset_A,Rset_b);
    
    Rvertex = Rset.V';
    
    % scaling factor constraint
    for i = 1:1:size(Rvertex,2)
        AEMPCineq = [AEMPCineq; -Omega.H(:,1:nx) zeros(size(Omega.H(:,1:nx),1),nx*(t+N)+nu*N) -Omega.H(:,nx+1) zeros(size(Omega.H(:,1:nx),1),N)];
        bEMPCineq = [bEMPCineq; -Omega.H(:,1:nx)*Rvertex(:,i)];
    end
    
%     AEMPCeq = [];
%     bEMPCeq = [];
%     for i = 1:N
%         AEMPCeq = [AEMPCeq; zeros(1,nx*(t+N+1)+N) zeros(1,i) 1 zeros(1,N-i)];
%         bEMPCeq = [bEMPCeq; 1];
%     end
    
    % output constraint 
    COmega = 0.8607;
    AEMPCineq = [AEMPCineq; [-C zeros(ny,nx*(t+N)) zeros(ny,N) -COmega zeros(ny,N)]; [C zeros(ny,nx*(t+N)) zeros(ny,N) -COmega zeros(ny,N)]];
    bEMPCineq = [bEMPCineq; -Y(t+1)+etah; Y(t+1)-etal];
    
    
    % optimization
    s0 = [kron(ones(t+N+1,1),Xesti(:,1)); kron(ones(N,1),0); kron(ones(N+1,1),1)];
%     s0 = [kron(ones(t+N+1,1),[0;0]); kron(ones(N,1),0); kron(ones(N+1,1),1)];
    
%     tic;
%     [s_vec, fval] = fmincon(@empc_obj,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);
%     t_fmincon(t+1) = toc;
%     flops_fmincon(t+1) = t_fmincon(t+1)/t_flop;

    problem = createOptimProblem('fmincon','objective',@empc_obj,'x0',s0,'Aineq',AEMPCineq,'bineq',bEMPCineq,'Aeq',AEMPCeq,'beq',bEMPCeq,'lb',lb,'ub',ub,'nonlcon',@empc_nonlcon,'options',options);
    gs = GlobalSearch;
    [s_vec, fval] = run(gs,problem);

    Z(:,1) = s_vec(1:nx);
    V(1) = s_vec(nx*(1+t+N)+1);
    U(1) = V(1) + K*(Xesti(:,1) - Z(:,1));
    Gamma(1) = s_vec(nx*(1+t+N)+1+N);
    X(:,2) = A*X(:,1) + B*U(1) + Bw*W(:,1);
    Xesti(:,2) = A*Xesti(:,1) + B*U(1) + L*(Y(1)-C*Xesti(:,1));
    Y(2) = C*X(:,2) + Eta(2);
    VavSum(1) = 50*X(1,1);
    AsyAve(1) = sum(VavSum)/1;
%     VavSumUpper(1) = Z(1,1)+Gamma(1)*0.018859+0.01*(Z(2,1)+Gamma(1)*0.86064795)^2;
    VavSumUpper(1) = 50*Z(1,1);
    for j = 1:size(Omegavertex,2)
        VavSumUpper(1) = max(VavSumUpper(1),50*(Z(1,1)+Gamma(1)*Omegavertex(1,j)));
    end
    AsyAveUpper(1) = sum(VavSumUpper)/1;
%     VavSumLower(1) = Z(1,1)-Gamma(1)*0.02704;
%     AsyAveLower(1) = sum(VavSumLower)/1;
%     
%     figure
%     plot(Z(:,t+1)+Gamma(t+1)*Omega);
%     hold on
%     Rset.plot('color','g')
%     hold off
    
    % t >= 1
    for t = 1:M-1
        t
%         s0 = [s_vec(1:nx*(t+N)); A*s_vec(nx*(t+N)-1:nx*(t+N))+B*(vs+K*(s_vec(nx*(t+N)-1:nx*(t+N))-zs)); s_vec(2+nx*(t+N):nx*(t+N)+N); vs+K*(s_vec(nx*(t+N)-1:nx*(t+N))-zs); s_vec(nx*(t+N)+1+N+1:end); 1];
        s0 = [s_vec(1:nx*(t+N)); zs; s_vec(2+nx*(t+N):nx*(t+N)+N); vs; s_vec(nx*(t+N)+1+N+1:end); 1];
        % bounds
        lb = [kron(ones(t+N+1,1),zl); kron(ones(N,1),vl); kron(ones(N+1,1),0)];
        ub = [kron(ones(t+N+1,1),zh); kron(ones(N,1),vh); kron(ones(N+1,1),1)];
        
        % initial state constraint
        AEMPCineq = [zeros(size(Xi.H(:,1:nx),1),t*nx) -Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N+N+1); zeros(size(Xi.H(:,1:nx),1),t*nx) Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N+N+1)];
        bEMPCineq = [1.15*Xi.H(:,nx+1) - Xi.H(:,1:nx)*Xesti(:,t+1); 1.15*Xi.H(:,nx+1) + Xi.H(:,1:nx)*Xesti(:,t+1)];
        
        % terminal equality constraint
        AEMPCeq = [zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N+N+1)];
        bEMPCeq = zs;
        
        % reachable set
        Rset = (Xesti(:,t+1)+E) & (A*Rset+B*U(t)+Bw*Pw);
        
        Rset_A = [Rset.H(:,1:nx); C; -C];
        Rset_b = [Rset.H(:,nx+1); Y(t+1)+etah; -(Y(t+1)+etal)];
        
        Rset = Polyhedron(Rset_A,Rset_b);
        Rvertex = Rset.V';
        
        if Rset.contains(X(:,t+1))==0
            error('contains');
        end
        
        % scaling factor constraint
        for i = 1:1:size(Rvertex,2)
            AEMPCineq = [AEMPCineq; zeros(size(Omega.H(:,1:nx),1),t*nx) -Omega.H(:,1:nx) zeros(size(Omega.H(:,1:nx),1),nx*N+nu*N) -Omega.H(:,nx+1) zeros(size(Omega.H(:,1:nx),1),N); zeros(size(Omega.H(:,1:nx),1),t*nx) Omega.H(:,1:nx) zeros(size(Omega.H(:,1:nx),1),nx*N+nu*N) -Omega.H(:,nx+1) zeros(size(Omega.H(:,1:nx),1),N)];
            bEMPCineq = [bEMPCineq; -Omega.H(:,1:nx)*Rvertex(:,i); Omega.H(:,1:nx)*Rvertex(:,i)];
        end
        
%         AEMPCeq = [];
%         bEMPCeq = [];
%         for i = 1:N
%             AEMPCeq = [AEMPCeq; zeros(1,nx*(t+N+1)+N) zeros(1,i) 1 zeros(1,N-i)];
%             bEMPCeq = [bEMPCeq; 1];
%         end
        
        % output constraint 
        for i = 0:1:t-1
            AEMPCineq = [AEMPCineq; [zeros(ny,i*nx) -C zeros(ny,nx*(t-i+N)) zeros(ny,N) zeros(ny,N+1)]; [zeros(ny,i*nx) C zeros(ny,nx*(t-i+N)) zeros(ny,N) zeros(ny,N+1)]];
            bEMPCineq = [bEMPCineq; Gamma(i+1)*COmega-Y(i+1)+etah; Gamma(i+1)*COmega+Y(i+1)-etal];
        end
        
        AEMPCineq = [AEMPCineq; [zeros(ny,t*nx) -C zeros(ny,nx*(N)) zeros(ny,N) -COmega zeros(ny,N)]; [zeros(ny,t*nx) C zeros(ny,nx*(N)) zeros(ny,N) -COmega zeros(ny,N)]];
        bEMPCineq = [bEMPCineq; -Y(t+1)+etah; Y(t+1)-etal];
                
        % optimization
        tic;
        [s_vec, fval, exitflag] = fmincon(@empc_obj,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);
        t_fmincon(t+1) = toc;
        flops_fmincon(t+1) = t_fmincon(t+1)/t_flop;
%         problem = createOptimProblem('fmincon','objective',@empc_obj,'x0',s0,'Aineq',AEMPCineq,'bineq',bEMPCineq,'Aeq',AEMPCeq,'beq',bEMPCeq,'lb',lb,'ub',ub,'nonlcon',@empc_nonlcon,'options',options);
%         gs = GlobalSearch;
%         [s_vec, fval, exitflag] = run(gs,problem);

        if exitflag<0
            problem = createOptimProblem('fmincon','objective',@empc_obj,'x0',s0,'Aineq',AEMPCineq,'bineq',bEMPCineq,'Aeq',AEMPCeq,'beq',bEMPCeq,'lb',lb,'ub',ub,'nonlcon',@empc_nonlcon,'options',options);
            gs = GlobalSearch;
            [s_vec, fval, exitflag] = run(gs,problem);
            if exitflag<0
                error('error')
            end
        end
        
        Z(:,t+1) = s_vec(t*nx+1:nx*(t+1));
        V(t+1) = s_vec(1+nx*(t+N+1));
        U(t+1) = V(t+1) + K*(Xesti(:,t+1) - Z(:,t+1));
        Gamma(t+1) = s_vec(nx*(1+t+N)+1+N);
        X(:,t+2) = A*X(:,t+1) + B*U(t+1) + Bw*W(:,t+1);
        Xesti(:,t+2) = A*Xesti(:,t+1) + B*U(t+1) + L*(Y(t+1)-C*Xesti(:,t+1));
        Y(t+2) = C*X(:,t+2) + Eta(t+2);
        VavSum(t+1) = 50*X(1,t+1);
        AsyAve(t+1) = sum(VavSum)/(t+1);
        VavSumUpper(t+1) = 50*Z(1,t+1);
        for j = 1:size(Omegavertex,2)
            VavSumUpper(t+1) = max(VavSumUpper(t+1),50*(Z(1,t+1)+Gamma(t+1)*Omegavertex(1,j)));
        end
        AsyAveUpper(t+1) = sum(VavSumUpper)/(t+1);
%         VavSumUpper(t+1) = Z(1,t+1)+Gamma(t+1)*0.018859+0.01*(Z(2,1)+Gamma(t+1)*0.86064795)^2;
%         AsyAveUpper(t+1) = sum(VavSumUpper)/(t+1);
%         VavSumLower(t+1) = Z(1,t+1)-Gamma(t+1)*0.01894;
%         AsyAveLower(t+1) = sum(VavSumLower)/(t+1);
        
%         figure
%         plot(Z(:,t+1)+Gamma(t+1)*Omega);
%         hold on
%         Rset.plot('color','g')
%         
%         figure
%         plot(Z(:,t+1)+Gamma(t+1)*Omega);
%         hold on;
%         plot(Z(:,t+1)+Xi);
%         plot(X(1,t+1),X(2,t+1),'k*');
%         plot(Xesti(1,t+1),Xesti(2,t+1),'bo');
        
    end
    
    
    Ax = [eye(nx); -eye(nx)];
    bx = [xh; -xl];
    Px = Polyhedron(Ax,bx);
    
    Az = [eye(nx); -eye(nx)];
    bz = [zh; -zl];
    Pz = Polyhedron(Az,bz);
    
    
    figure
    s0 = Px.plot('color','r');
    hold on
    s1 = Pz.plot('color','g');
    
    for i = 1:1:M
        set_omega=Z(:,i)+Gamma(i)*Omega;
        set_omega.plot('alpha',0.3,'color',[0.3 0.3 0.3]);
    end
    
    
    p1 = plot(X(1,:),X(2,:),'LineWidth',1.5);
    p2 = plot(Xesti(1,:),Xesti(2,:),'--','LineWidth',1.5);
    p3 = plot(Z(1,:),Z(2,:),':','LineWidth',1.5);
    
    grid on
    xlabel({'$x_{1}$'},'Interpreter','latex');
    ylabel({'$x_{2}$'},'Interpreter','latex');
    legend([s0 s1 p1, p2, p3],{'$X$','$Z$','$x$','$\hat{x}$','$z$'},'Interpreter','latex','Location','bestoutside');
    
    figure
    subplot(2,1,1)
    plot(0:M-1,U,0:M-1,V,'--','LineWidth',1.5);
    grid on;
    legend({'$u$','$v$'},'Interpreter','Latex');
    legend boxoff
    xlim([0 M])
%     ylim([-1 0.5])
    subplot(2,1,2)
    plot(0:M-1,Gamma,'LineWidth',1.5);
    hold on
    plot(0:M-1,Gamma,'--','LineWidth',1.5);
    grid on;
    legend({'$\gamma$','$\underline{\gamma}$'},'Interpreter','Latex');
    legend boxoff
    xlim([0 M])
%     ylim([0 1.2]);
    set(gca,'YTick',[0 0.5 1 1.2]);
    
    
    
    
    
    
    
    fprintf('average time for fmincon: %g\tflops = %g\n', mean(t_fmincon), mean(flops_fmincon));

    
    
%     %%%%  sampling time = 0.1 min
%     Ts = 0.1;
%     figure
%     subplot(2,1,1)
%     plot(0:Ts:(M-1)*Ts,U,0:Ts:(M-1)*Ts,V,'--','LineWidth',1.5);
%     grid on;
%     legend({'$u$','$v$'},'Interpreter','Latex');
%     legend boxoff
%     ylabel('Manipulated variable')
%     subplot(2,1,2)
%     plot(0:Ts:(M-1)*Ts,Gamma,'LineWidth',1.5);
%     hold on
%     plot(0:Ts:(M-1)*Ts,Gamma,'--','LineWidth',1.5);
%     grid on;
%     legend({'$\gamma$','$\underline{\gamma}$'},'Interpreter','Latex');
%     legend boxoff
%     ylabel('Scaling factor')
%     xlabel('Time (min)')
%     ylim([0 1.2]);
    
    
    
    
    
    