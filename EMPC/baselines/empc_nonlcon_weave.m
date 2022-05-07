function [c,ceq] = empc_nonlcon_weave(s)
    global t N V A B;
    nx = 2;
    ceq = [];
    %% load disturbance
    load("multidisturbances.mat");

    %% for the first scenario
    Xesti1 = s(1:nx*(t+1)); % state at current timestep
    Xpred1 = s(nx*t+1:nx*(t+N+1)); % states for future timesteps (N+1)
    Vvec1 = s(nx*(1+t+N)+1:nx*(t+N+1)+N); % inputs for future timesteps (N)
    % prediction
    for i = 1:1:N
        ceq = [ % in the end would have size of 24x1 -- 2*N
            ceq;
            A*Xpred1(nx*i-1:nx*i)+B*Vvec1(i)+d1(:,t+i) - Xpred1(nx*(i+1)-1:nx*(i+1)) % the difference in states between future timesteps and prediction
            ];
    end
    % estimation
    if t>0
        for i = 1:1:t
            ceq = [ceq; A*Xesti1(nx*i-1:nx*i)+B*V(i) - Xesti1(nx*(i+1)-1:nx*(i+1))];
        end
    end

    %% for the second scenario
    bias = nx*(t+N+1)+N;
    Xesti2 = s(bias+1:bias+nx*(t+1));
    Xpred2 = s(bias+nx*t+1:bias+nx*(t+N+1));
    Vvec2 = s(bias+nx*(1+t+N)+1:bias+nx*(t+N+1)+N);
    % prediction
    for i = 1:1:N
        ceq = [ % in the end would have size of 24x1 -- 2*N
            ceq;
            A*Xpred2(nx*i-1:nx*i)+B*Vvec2(i)+d2(:,t+i) - Xpred2(nx*(i+1)-1:nx*(i+1)) % the difference in states between future timesteps and prediction
            ];
    end
    % estimation
    if t>0
        for i = 1:1:t
            ceq = [ceq; A*Xesti2(nx*i-1:nx*i)+B*V(i) - Xesti2(nx*(i+1)-1:nx*(i+1))];
        end
    end
    
    %% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    c = [];
end
