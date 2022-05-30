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
% d2 is the nominal case
d1(:,:) = 0.08;
d2(:,:) = 0;


save('multidisturbances.mat','d1','d2','d3');