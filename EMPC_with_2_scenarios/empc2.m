clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    EMPC algorithm    %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    format long

    global t N A B zs Omegavertex w1 w2;
    w1 = 0.15;
    w2 = 0.85;
    
    % system description
    A = [1 1; 0 1];
    B = [0.5; 1];
    C = [0 1];
    K = [-0.6696 -1.3261];
    
    load('Xi.mat')
    Omegavertex = Xi.V';
    
    nx = 2;
    nu = 1;
    ny = 1;
    
    % constraints
    xh = [8; 10];
    xl = [-8; -10];
    
    uh = 10;
    ul = -10;
    
    wh = [1; 1];
    wl = [-1; -1];
    
    
    Aw = [eye(2); -eye(2)];
    bw = [wh; -wl];
    Pw = Polyhedron(Aw,bw);
    
    % steady state   
    zs = [0;0];
    vs = 0;
    % the tightened constraints
    % for states
    zh = [5.48373; 7.5];
    zl = [-5.48373; -7.5];
    
    % for inputs
    vh = 5;
    vl = -5;
    
    % horizons
    M = 50;
    N = 12;
    
    
    % state disturbance
    load('disturbances.mat');

    % input and output data
    V = zeros(nu,M);
    U = zeros(nu,M);
    Y = zeros(ny,M+1);
    VavSum = zeros(M,1);
    AsyAve = zeros(M,1);
    
    % initial state
    X = [[-0.2; 2.5] zeros(nx,M)];
    Z = zeros(nx,M);
    Y(1) = C*X(:,1);
    
    
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10);

    t = 0;
    
    % bound
    lb = [
        kron(ones(N+1,1),zl);
        kron(ones(N,1),vl);
        kron(ones(N+1,1),zl);
        kron(ones(N,1),vl)
        ];
    ub = [
        kron(ones(N+1,1),zh);
        kron(ones(N,1),vh);
        kron(ones(N+1,1),zh);
        kron(ones(N,1),vh)
        ];
    
    % initial state constraint
    AEMPCineq = [
        -Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N) zeros(size(-Xi.H(:,1:nx))) zeros(size(Xi.H,1),nx*(t+N)+N);
        Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N) zeros(size(Xi.H(:,1:nx))) zeros(size(Xi.H,1),nx*(t+N)+N);
        zeros(size(-Xi.H(:,1:nx))) zeros(size(Xi.H,1),nx*(t+N)+N) -Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N);
        zeros(size(Xi.H(:,1:nx))) zeros(size(Xi.H,1),nx*(t+N)+N) Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N)
        ];
    bEMPCineq = [
        Xi.H(:,nx+1)-Xi.H(:,1:nx)*X(:,t+1);
        Xi.H(:,nx+1)+Xi.H(:,1:nx)*X(:,t+1);
        Xi.H(:,nx+1)-Xi.H(:,1:nx)*X(:,t+1);
        Xi.H(:,nx+1)+Xi.H(:,1:nx)*X(:,t+1)
        ];
    
    
    % terminal equality constraint
    AEMPCeq = [
        zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N) zeros(nx,nx*(t+N)) zeros(size(eye(nx))) zeros(nx,N);
        zeros(nx,nx*(t+N)) zeros(size(eye(nx))) zeros(nx,N) zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N);
        eye(nx) zeros(nx,nx*(t+N)) zeros(nx,N) zeros(nx,nx*(t+N+1)) zeros(nx,N);
        zeros(nx,nx*(t+N+1)) zeros(nx,N) eye(nx) zeros(nx,nx*(t+N)) zeros(nx,N);
        zeros(1,nx*(t+N+1)) 1 zeros(1,N-1) zeros(1,nx*(t+N+1)) -1 zeros(1,N-1)
        ];
    bEMPCeq = [
        zs;
        zs;
        X(:,1);
        X(:,1);
        0
        ];
    
    
    
    s0 = [
        kron(ones(N+1,1),X(:,1));
        kron(ones(N,1),0);
        kron(ones(N+1,1),X(:,1));
        kron(ones(N,1),0)
        ];
    
    [s_vec, fval] = fmincon(@empc_obj2,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon2,options);

    Z(:,1) = s_vec(1:nx);
    V(1) = s_vec(nx*(1+t+N)+1);
    U(1) = V(1) + K*(X(:,1) - Z(:,1));
    X(:,2) = A*X(:,1) + B*U(1) + W(:,1);
    Y(2) = C*X(:,2);
    VavSum(1) = X(1,1)^2+X(2,1)^2;
    AsyAve(1) = sum(VavSum)/1;
    
    % t >= 1
    for t = 1:M-1
        t
        s0 = [s_vec(nx+1:nx*(1+N));
            A*s_vec(nx*N+1:nx*(1+N))+B*(vs+K*(s_vec(nx*N+1:nx*(1+N))-zs));
            s_vec(2+nx*(1+N):nx*(1+N)+N);
            vs+K*(s_vec(nx*N+1:nx*(1+N))-zs);
            s_vec(nx+1:nx*(1+N));
            A*s_vec(nx*N+1:nx*(1+N))+B*(vs+K*(s_vec(nx*N+1:nx*(1+N))-zs));
            s_vec(2+nx*(1+N):nx*(1+N)+N);
            vs+K*(s_vec(nx*N+1:nx*(1+N))-zs)
            ];

        
        % initial state constraint
        AEMPCineq = [
            -Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N) zeros(size(-Xi.H(:,1:nx))) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N);
            Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N) zeros(size(Xi.H(:,1:nx))) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N);
            zeros(size(-Xi.H(:,1:nx))) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N) -Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N);
            zeros(size(Xi.H(:,1:nx))) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N) Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N)
            ];
        bEMPCineq = [
            Xi.H(:,nx+1) - Xi.H(:,1:nx)*X(:,t+1);
            Xi.H(:,nx+1) + Xi.H(:,1:nx)*X(:,t+1);
            Xi.H(:,nx+1) - Xi.H(:,1:nx)*X(:,t+1);
            Xi.H(:,nx+1) + Xi.H(:,1:nx)*X(:,t+1)
            ];
        
        % terminal equality constraint
        AEMPCeq = [
            zeros(nx,nx*N) eye(nx) zeros(nx,N) zeros(nx,nx*N) zeros(size(eye(nx))) zeros(nx,N);
            zeros(nx,nx*N) zeros(size(eye(nx))) zeros(nx,N) zeros(nx,nx*N) eye(nx) zeros(nx,N);
            eye(nx) zeros(nx,nx*N) zeros(nx,N) zeros(nx,nx*(N+1)) zeros(nx,N);
            zeros(nx,nx*(N+1)) zeros(nx,N) eye(nx) zeros(nx,nx*N) zeros(nx,N);
            zeros(1,nx*(N+1)) 1 zeros(1,N-1) zeros(1,nx*(N+1)) -1 zeros(1,N-1)
            ];
        bEMPCeq = [
            zs;
            zs;
            X(:,t+1);
            X(:,t+1);
            0
            ];
                
        % optimization
        [s_vec, fval, exitflag] = fmincon(@empc_obj2,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon2,options);

%         if exitflag<0
%             problem = createOptimProblem('fmincon','objective',@empc_obj2,'x0',s0,'Aineq',AEMPCineq,'bineq',bEMPCineq,'Aeq',AEMPCeq,'beq',bEMPCeq,'lb',lb,'ub',ub,'nonlcon',@empc_nonlcon2,'options',options);
%             gs = GlobalSearch;
%             [s_vec, fval, exitflag] = run(gs,problem);
%             if exitflag<0
%                 error('error')
%             end
%         end
        
        Z(:,t+1) = s_vec(1:nx);
        V(t+1) = s_vec(1+nx*(N+1));
        U(t+1) = V(t+1) + K*(X(:,t+1) - Z(:,t+1));
        X(:,t+2) = A*X(:,t+1) + B*U(t+1) + W(:,t+1);
        Y(t+2) = C*X(:,t+2);
        VavSum(t+1) = X(1,t+1)^2+X(2,t+1)^2;
        AsyAve(t+1) = sum(VavSum)/(t+1);  
    end
    
    
    Ax = [eye(nx); -eye(nx)];
    bx = [xh; -xl];
    Px = Polyhedron(Ax,bx);
    
    Az = [eye(nx); -eye(nx)];
    bz = [zh; -zl];
    Pz = Polyhedron(Az,bz);
    
    
%     figure
%     s0 = Px.plot('color','r');
%     hold on
%     s1 = Pz.plot('color','g');
    
%     for i = 1:1:M
%         set_omega=Z(:,i)+Xi;
%         set_omega.plot('alpha',0.3,'color',[0.3 0.3 0.3]);
%     end
    
    
%     p1 = plot(X(1,:),X(2,:),'LineWidth',1.5);
%     p3 = plot(Z(1,:),Z(2,:),':','LineWidth',1.5);
    
%     grid on
%     xlabel({'$x_{1}$'},'Interpreter','latex');
%     ylabel({'$x_{2}$'},'Interpreter','latex');
%     legend([s0 s1 p1, p3],{'$X$','$Z$','$x$','$z$'},'Interpreter','latex','Location','bestoutside');
    
    figure;
    plot(1:M,VavSum);
    hold on;
    plot(1:M,AsyAve);
    legend('cost','average cost')
