clear;
clc;
% The disturbance: W = {ω ∈ R2 : |ω1| ≤ 0.2, |ω2| ≤ 0.05}
% |Eta| <= 0.01
% time steps T = 100;
T = 50;

%% gennerate the noise
%% uniformlly distributed within the bound
% w1 = (rand(1,T)-0.5)/2.5;
% w2 = (rand(1,T)-0.5)/10;
% Eta = (rand(1,T+1)-0.5)/50;
% W = [w1;w2];

%% Noise are all zero
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% Eta = zeros(1,T+1);
% W = [w1;w2];

%% normally distributed noise centered at zero bounded within the restriction
% p = 6;
% Eta = -0.01 + (0.01+0.01)*sum(rand(T+1,p),2)/p;
% w1 = -0.2 + (0.2+0.2)*sum(rand(T,p),2)/p;
% w2 = -0.05 + (0.05+0.05)*sum(rand(T,p),2)/p;
% W = [w1';w2'];


%% with large noise in the beginning(first 25 samples)
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% w1(1:25) = 0.2;
% w2(1:25) = 0.05;
% Eta = zeros(1,T+1);
% Eta(1:25) = 0.01;
% W = [w1;w2];

%% with large noise in middle(45:70)
w1 = zeros(1,T);
w2 = zeros(1,T);
w1(12:25) = 1;
w2(12:25) = -0.3;
w1(30:32) = -0.8;
w2(30:32) = 0.02;
w1(1:5) = 0.5;
w2(1:5) = 0.2;
W = [w1;w2];

%% with large noise in the end(75:100,75:101)
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% w1(75:T) = 0.1;
% w2(75:T) = 0.1;
% Eta = zeros(1,T+1);
% Eta(75:T+1) = 0.01;
% W = [w1;w2];



save('disturbances.mat','W');


%     if(t >= 2 && VavSum(t) > VavSum(t-1))
%         % if the current cost is larger than the previous
%         % we need to modify our input U at time t+1
%         if(sum(X(:,t+1) > sum(X(:,t))))
%             % if the current X is larger, for this cost function, we
%             % need to reduce X to reduce the cost
%             % This can be done by simply substract a bias from X
%             bias = [0.01;0.4];
%             X(:,t+2) = A*X(:,t+1) + B*U(t+1) + Bw*W(:,t+1) - bias;
%             Xesti(:,t+2) = A*Xesti(:,t+1) + B*U(t+1) + L*(Y(t+1)-C*Xesti(:,t+1)) - bias;
%             Y(t+2) = C*X(:,t+2) + Eta(t+2);
%         end
%     end