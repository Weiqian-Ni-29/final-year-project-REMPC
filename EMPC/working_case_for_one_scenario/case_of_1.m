clear;
close all;
clc;

%% Rigid tube based EMPC algorithm, case of one disturbance

global t N V A B zs w;
w = 0.8;

%% System discription
A = [1 1; 0 1];
B = [0.5; 1];
C = [0 1];
K = [-0.6696 -1.3261];

% load the RPI set
load('Xi.mat')

% size of states, input and output
nx = 2;
nu = 1;
ny = 1;

%% constraints
% for real states 
xh = [8; 5];
xl = [-8; -5];

% for real inputs
uh = 3;
ul = -3;

% the bound for disturbances
wh = [0.1; 0.1];
wl = [-0.1; -0.1];

% load the steady state of operation
load('newSteadyState.mat');

% the tightened constraints
% for states
zh = [7.5732; 4.725];
zl = [-7.5732; -4.725];

% for inputs
vh = 2.34956;
vl = -2.34956;


%% the horizons
T = 50; % the overall simulation timestep
N = 12; % the prediction horizon

%% load the disturbances
% the real disturbances used in the simulation
load('disturbances.mat');

% the predicted disturbances
load('multidisturbances.mat');

%% initialize the input and output, costs and states
V = zeros(nu,T);
U = zeros(nu,T);
Y = zeros(ny,T+1);

stage_cost = zeros(T,1);
average_cost = zeros(T,1);

X = [[0.2; -3] zeros(nx,T)];
Z = zeros(nx,T);
Y(1) = C * X(:,1);

% the bound for state and input;
lb = [kron(ones(N+1,1),zl); kron(ones(N,1),vl)];
ub = [kron(ones(N+1,1),zh); kron(ones(N,1),vh)];

% specify the optimoption
options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-10,'TolFun',1e-10,'TolX',1e-10);

% initialize the time t = 0;
t = 0;

% set the initialize state constraint
% force the initial state be within the RPI set
AEMPCineq = [-Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*N+N); Xi.H(:,1:nx) zeros(size(Xi.H,1),nx*N+N)];
bEMPCineq = [1.03*Xi.H(:,nx+1)-Xi.H(:,1:nx)*X(:,t+1); 1.03*Xi.H(:,nx+1)+Xi.H(:,1:nx)*X(:,t+1)];

% output constraint
COmega = 0.4268;
AEMPCineq = [
    AEMPCineq;
    [-C zeros(ny,nx*N) zeros(ny,N)];
    [C zeros(ny,nx*N) zeros(ny,N)]
    ];
bEMPCineq = [
    bEMPCineq;
    -Y(t+1)+COmega;
    Y(t+1)+COmega
];

% make the constraints for the disturbance
disturbance_constraint = [];
for i = 0:N-1
    disturbance_constraint = [
        disturbance_constraint;
        [zeros(nx,i) A -eye(nx) zeros(nx,nx*(N-1)-i) zeros(nx,i) B zeros(nx,N-1-i)]
        ];
end

disturbances = [];
for i = 1:N
    disturbances = [
        disturbances;
        d1(:,i);
    ];
end

% equality constraint
AEMPCeq = [
    zeros(nx,nx*N) eye(nx) zeros(nx,N); % terminal equality constraint
    eye(nx) zeros(nx,nx*N) zeros(nx,N); % the first state constraint
    %disturbance_constraint % constraints for disturbance
    ];
bEMPCeq = [
    zs;
    [0.2;-3];
    %disturbances * (-1)
    ]; % terminal equality constraint;

% initialize the first state and inputs
s0 = [
    kron(ones(N+1,1),X(:,1)); % the first and next N states
    kron(ones(N,1),0) %  the N inputs for the N predicted states
];
% the predicted disturbances are implemented in the nonlinear constraints
[s_vec, fval] = fmincon(@empc_obj,s0,[],[],AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);

Z(:,1) = s_vec(1:nx);
V(1) = s_vec(nx*(1+N)+1);
U(1) = V(1) + K*(X(:,1) - Z(:,1));
X(:,2) = A*X(:,1) + B*U(1) + W(:,1);
Y(2) = C*X(:,2);
stage_cost(1) = X(1,1)+0.01*X(2,1)^2;
average_cost(1) = sum(stage_cost)/1;

for t = 1:T-1
    t

    % make the constraints for the disturbance
    disturbances = [];
    for i = t+1:N+t
        disturbances = [
            disturbances;
            d1(:,i);
        ];
    end

    % equality constraint
    AEMPCeq = [
        zeros(nx,nx*N) eye(nx) zeros(nx,N); % terminal equality constraint
        eye(nx) zeros(nx,nx*N) zeros(nx,N); % ensures the nominal state is moving to the next step
        %disturbance_constraint
        ];
    bEMPCeq = [
        zs; % terminal equality constraint;
        s_vec(3); % ensures the nominal state is moving to the next step
        s_vec(4);
        %disturbances*(-1)
        ]; 
    s0 = [
        s_vec(3:end);
        0;
        0
        ];
    % the bounds and constraints are all the same for the base case
    [s_vec, fval] = fmincon(@empc_obj,s0,[],[],AEMPCeq,bEMPCeq,lb,ub,@empc_nonlcon,options);
    Z(:,t+1) = s_vec(1:nx);
    V(t+1) = s_vec(nx*(1+N)+1);
    U(t+1) = V(t+1) + K*(X(:,t+1) - Z(:,t+1));
    X(:,t+2) = A*X(:,t+1) + B*U(t+1) + W(:,t+1);
    Y(t+2) = C*X(:,t+2);
    stage_cost(t+1) = X(1,t+1)+0.01*X(2,t+1)^2;
    average_cost(t+1) = sum(stage_cost)/(t+1);
end

%% plots
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

for i = 1:1:T
    set_omega=Z(:,i)+Xi;
    set_omega.plot('alpha',0.3,'color',[0.3 0.3 0.3]);
end

p1 = plot(X(1,:),X(2,:),'LineWidth',1.5);
p3 = plot(Z(1,:),Z(2,:),':','LineWidth',1.5);

grid on;
xlabel({'$x_{1}$'},'Interpreter','latex');
ylabel({'$x_{2}$'},'Interpreter','latex');
legend([s0 s1 p1, p3],{'$X$','$Z$','$x$','$z$'},'Interpreter','latex','Location','bestoutside');
figure;
plot(1:T,stage_cost);
hold on;
plot(1:T,average_cost);
legend('cost','average cost')
