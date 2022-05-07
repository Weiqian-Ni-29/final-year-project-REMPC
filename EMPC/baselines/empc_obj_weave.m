function [fval] = empc_obj_weave(s)
    global t N Omegavertex;
    nx = 2;
    %% load the two different disturbance scenarios
    load('multidisturbances.mat');
    %% declare the weights
    w1 = 0.5;
    w2 = 0.5;
    %% calculate the first cost function
    Xpred1 = s(nx*t+1:nx*(t+N+1));
    fval1 = 0;
    for i = 0:N-1
        f = (Xpred1(nx*i+1)+d1(1,t+i+1)) + 0.01*(Xpred1(nx*i+nx)+d1(2,t+i+1))^2;
        fval1 = fval1 + f;
    end
    
    %% calculate the second cost function
    Xpred2 = s(nx*(t+N+1)+N+1:nx*(t+N+1)+N+nx*(t+N+1)); 
    fval2 = 0;
    for i = 0:N-1
        f = (Xpred2(nx*i+1)+d2(1,t+i+1)) + 0.01*(Xpred2(nx*i+nx)+d2(2,t+i+1))^2;
        fval2 = fval2 + f;
    end
    fval = w1*fval1+w2*fval2;
end
