clear;
clc;
% The disturbance: W = {ω ∈ R2 : |ω1| ≤ 1, |ω2| ≤ 1}
% |Eta| <= 0.01
% time steps T = 100;
T = 50;

%% gennerate the noise
%% uniformlly distributed within the bound
% w1 = (rand(1,T)-0.5)/5;
% w2 = (rand(1,T)-0.5)/5;
% W = [w1;w2];

%% Noise are all zero
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% W = [w1;w2];

%% normally distributed noise centered at zero bounded within the restriction
% p = 6;
% w1 = -0.1 + (0.1+0.1)*sum(rand(T,p),2)/p;
% w2 = -0.1 + (0.1+0.1)*sum(rand(T,p),2)/p;
% W = [w1';w2'];


%% normally distributed noise centered at zero bounded within the restriction with bias
% p = 6;
% bias = 0.05;
% w1 = -0.1 + (0.1+0.1)*sum(rand(T,p),2)/p + bias;
% w2 = -0.1 + (0.1+0.1)*sum(rand(T,p),2)/p + bias;
% W = [w1';w2'];

%% with large noise in the beginning(first 25 samples)
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% w1(1:25) = 0.2;
% w2(1:25) = 0.05;
% W = [w1;w2];

%% with large noise in middle(45:70)
w1 = zeros(1,T);
w2 = zeros(1,T);
w1(12:25) = 1;
w2(12:25) = -0.3;
w1(30:32) = -0.8;
w2(30:32) = 0.02;
W = [w1;w2];

%% with large noise in the end(75:100,75:101)
% w1 = zeros(1,T);
% w2 = zeros(1,T);
% w1(75:T) = 0.1;
% w2(75:T) = 0.1;
% W = [w1;w2];



save('disturbances.mat','W');

