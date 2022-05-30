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
% d1 is the nominal case
d1(1,12:25) = 1;
d2(1,12:25) = -0.3;
d2(1,30:32) = -0.8;
d2(2,30:32) = 0.02;


save('multidisturbances.mat','d1','d2','d3');