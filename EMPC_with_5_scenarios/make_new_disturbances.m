clear;
clc;
T = 50;
global N;
%% scenario 1: no disturbance, scenario 2: no disturbance, scenario 3: no disturbance
% d1 = zeros(2,T+1+N);
% d2 = zeros(2,T+1+N);
% d3 = zeros(2,T+1+N);

%% scenario 1: a constant in the middle, scenario 3: no disturbance, scenario 2: a bump in the end
d1 = zeros(2,T+1+N);
d2 = zeros(2,T+1+N);
d3 = zeros(2,T+1+N);
d4 = zeros(2,T+1+N);
d5 = zeros(2,T+1+N);
% d1 is the nominal case
d2(1,:) = 1;
d2(2,:) = -0.3;
d3(1,:) = -0.8;
d3(2,:) = 0.02;
d4(1,:) = 0.5;
d4(2,:) = -0.2;
d5(1,:) = 1;
d5(2,:) = 1;

save('multidisturbances.mat','d1','d2','d3','d4','d5');