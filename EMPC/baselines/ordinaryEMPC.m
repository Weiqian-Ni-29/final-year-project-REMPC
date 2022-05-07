clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   Rigid Tube-based EMPC algorithm   %%%
%%% This is the most fundmental Tube-based
%%% EMPC algorithm, do not delete or change 
%%% this script(with the second model!!!), use this as a templete only!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

global t N V A B zs Omegavertex;

% system description
A = [1 1; 0 1];
B = [0.5;1];
C = [0 1];
K = [-0.6696 -1.3261];
load('Xi.mat');
Omegavertex = Xi.V';

nx = 2;
nu = 1;
ny = 1;

% constraints
xh = [8; 5];
xl = [-8; -5];

uh = 3;
ul = -3;

wh = [0.1; 0.1];
wl = [-0.1; -0.1];

Aw = [eye(2); -eye(2)];
bw = [wh; -wl];
Pw = Polyhedron(Aw,bw);

% steady state   
load('newSteadyState.mat')

% tightened constranits
zh = [7.5732; 4.725];
zl = [-7.5732; -4.725];

vh = 2.34956;
vl = -2.34956;

% horizons
M = 100;
N = 12;


% state disturbance
load('disturbances.mat')

% input and output data
V = zeros(nu,M);
U = zeros(nu,M);
Y = zeros(ny,M+1);
VavSum = zeros(M,1);
AsyAve = zeros(M,1);
VavSumUpper = zeros(M,1);
AsyAveUpper = zeros(M,1);

% initial state
X = [[0.2; -3] zeros(nx,M)];
Z = zeros(nx,M);
Y(1) = C*X(:,1);

options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10);

t = 0;

% bound
lb = [kron(ones(t+N+1,1),zl); kron(ones(N,1),vl)];
ub = [kron(ones(t+N+1,1),zh); kron(ones(N,1),vh)];

% initial state constraint
% the state cannot go outside of the RPI set
AEMPCineq = [-Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N); Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*(t+N)+N)];
bEMPCineq = [Xi.H(:,nx+1)-Xi.H(:,1:nx)*X(:,t+1); Xi.H(:,nx+1)+Xi.H(:,1:nx)*X(:,t+1)];


% terminal equality constraint
% force the last state (N+1) to be steady state
AEMPCeq = [zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N)];
bEMPCeq = zs;

% output constraint 
COmega = 0.4268;
% add additional two rows in the bottom
AEMPCineq = [AEMPCineq; [-C zeros(ny,nx*(t+N)) zeros(ny,N)]; [C zeros(ny,nx*(t+N)) zeros(ny,N)]];
bEMPCineq = [bEMPCineq; -Y(t+1)+COmega; Y(t+1)+COmega];


% optimization
s0 = [kron(ones(t+N+1,1),X(:,1)); kron(ones(N,1),0)];

[s_vec, fval] = fmincon(@empc_obj,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);

Z(:,1) = s_vec(1:nx);
V(1) = s_vec(nx*(1+t+N)+1);
U(1) = V(1) + K*(X(:,1) - Z(:,1));
X(:,2) = A*X(:,1) + B*U(1) + W(:,1);
Y(2) = C*X(:,2);
VavSum(1) = X(1,1)+0.01*X(2,1)^2;
AsyAve(1) = sum(VavSum)/1;

VavSumUpper(1) = Z(1,1) + 0.01*Z(2,1)^2;
for j = 1:size(Omegavertex,2)
    VavSumUpper(1) = max(VavSumUpper(1),Z(1,1)+Omegavertex(1,j) + 0.01*(Z(2,1)+Omegavertex(2,j))^2);
end
AsyAveUpper(1) = sum(VavSumUpper)/1;


% t >= 1
for t = 1:M-1
    t
    s0 = [s_vec(1:nx*(t+N)); A*s_vec(nx*(t+N)-1:nx*(t+N))+B*(vs+K*(s_vec(nx*(t+N)-1:nx*(t+N))-zs)); s_vec(2+nx*(t+N):nx*(t+N)+N); vs+K*(s_vec(nx*(t+N)-1:nx*(t+N))-zs)];

    % bounds
    lb = [kron(ones(t+N+1,1),zl); kron(ones(N,1),vl)];
    ub = [kron(ones(t+N+1,1),zh); kron(ones(N,1),vh)];
    
    % initial state constraint
    AEMPCineq = [zeros(size(Xi.H(:,1:nx),1),t*nx) -Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N); zeros(size(Xi.H(:,1:nx),1),t*nx) Xi.H(:,1:nx) zeros(size(Xi.H(:,1:nx),1),nx*N+nu*N)];
    bEMPCineq = [1.03*Xi.H(:,nx+1) - Xi.H(:,1:nx)*X(:,t+1); 1.03*Xi.H(:,nx+1) + Xi.H(:,1:nx)*X(:,t+1)];
    
    % terminal equality constraint
    AEMPCeq = [zeros(nx,nx*(t+N)) eye(nx) zeros(nx,N)];
    bEMPCeq = zs;
    
    % output constraint 
    for i = 0:1:t
        AEMPCineq = [AEMPCineq; [zeros(ny,i*nx) -C zeros(ny,nx*(t-i+N)) zeros(ny,N)]; [zeros(ny,i*nx) C zeros(ny,nx*(t-i+N)) zeros(ny,N)]];
        bEMPCineq = [bEMPCineq; COmega-Y(i+1); COmega+Y(i+1)];
    end
            
    % optimization

    [s_vec, fval, exitflag] = fmincon(@empc_obj,s0,AEMPCineq,bEMPCineq,AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);

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
    U(t+1) = V(t+1) + K*(X(:,t+1) - Z(:,t+1));
    X(:,t+2) = A*X(:,t+1) + B*U(t+1) + W(:,t+1);
    Y(t+2) = C*X(:,t+2);
    VavSum(t+1) = X(1,t+1)+0.01*X(2,t+1)^2;
    AsyAve(t+1) = sum(VavSum)/(t+1);
    VavSumUpper(t+1) = Z(1,t+1) + 0.01*Z(2,t+1)^2;
    for j = 1:size(Omegavertex,2)
        VavSumUpper(t+1) = max(VavSumUpper(t+1),Z(1,t+1)+Omegavertex(1,j)+0.01*(Z(2,t+1)+Omegavertex(2,j))^2);
    end
    AsyAveUpper(t+1) = sum(VavSumUpper)/(t+1);  
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
    set_omega=Z(:,i)+Xi;
    set_omega.plot('alpha',0.3,'color',[0.3 0.3 0.3]);
end


p1 = plot(X(1,:),X(2,:),'LineWidth',1.5);
p3 = plot(Z(1,:),Z(2,:),':','LineWidth',1.5);

grid on
xlabel({'$x_{1}$'},'Interpreter','latex');
ylabel({'$x_{2}$'},'Interpreter','latex');
legend([s0 s1 p1, p3],{'$X$','$Z$','$x$','$z$'},'Interpreter','latex','Location','bestoutside');
figure;
plot(1:100,VavSum);
hold on;
plot(1:100,AsyAve);
legend('cost','average cost')
