function [c,ceq] = empc_nonlcon2(s)
    global t N V A B w1 w2;
    
    nx = 2;
    Xpred1 = s(1:nx*(N+1)); % states for future timesteps (N+1), scenario 1
    Vvec1 = s(nx*(1+N)+1:nx*(N+1)+N); % inputs for future timesteps (N), scenario 1
    Xpred2 = s(nx*(N+1)+N+1:nx*(N+1)+N+nx*(N+1)); % states for future timesteps (N+1), scenario 2
    Vvec2 = s(nx*(N+1)+N+nx*(N+1)+1:nx*(N+1)+N+nx*(N+1)+N); % inputs for future timesteps (N), scenario 2
    load('multidisturbances.mat')
    ceq = [];
    % prediction
    %% for scenario 1
    for i = 1:1:N
        ceq = [ 
            ceq;
            A*Xpred1(nx*i-1:nx*i)+B*Vvec1(i) + d1(:,t+i) - Xpred1(nx*(i+1)-1:nx*(i+1)) % the difference in states between future timesteps and prediction
            ];
    end

    %% for scenario 2
    for i = 1:1:N
        ceq = [
            ceq;
            A*Xpred2(nx*i-1:nx*i)+B*Vvec2(i) + d2(:,t+i) - Xpred2(nx*(i+1)-1:nx*(i+1))
        ];
    end

    %% output 
    c = [];
end