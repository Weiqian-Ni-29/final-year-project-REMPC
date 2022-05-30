function [c,ceq] = empc_nonlcon5(s)
    global t N V A B;
    
    nx = 2;
    Xpred1 = s(1:nx*(N+1)); % states for future timesteps (N+1), scenario 1
    Vvec1 = s(nx*(1+N)+1:nx*(N+1)+N); % inputs for future timesteps (N), scenario 1
    Xpred2 = s(nx*(N+1)+N+1:nx*(N+1)+N+nx*(N+1)); % states for future timesteps (N+1), scenario 2
    Vvec2 = s(nx*(N+1)+N+nx*(N+1)+1:nx*(N+1)+N+nx*(N+1)+N); % inputs for future timesteps (N), scenario 2
    Xpred3 = s(2*(nx*(N+1)+N)+1:2*(nx*(N+1)+N)+nx*(N+1));
    Vvec3 = s(2*(nx*(N+1)+N)+nx*(N+1)+1:2*(nx*(N+1)+N)+nx*(N+1)+N);
    Xpred4 = s(3*(nx*(N+1)+N)+1:3*(nx*(N+1)+N)+nx*(N+1));
    Vvec4 = s(3*(nx*(N+1)+N)+nx*(N+1)+1:3*(nx*(N+1)+N)+nx*(N+1)+N);
    Xpred5 = s(4*(nx*(N+1)+N)+1:4*(nx*(N+1)+N)+nx*(N+1));
    Vvec5 = s(4*(nx*(N+1)+N)+nx*(N+1)+1:4*(nx*(N+1)+N)+nx*(N+1)+N);

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

    %% for scenario 3
    for i = 1:1:N
        ceq = [
            ceq;
            A*Xpred3(nx*i-1:nx*i)+B*Vvec3(i) + d3(:,t+i) - Xpred3(nx*(i+1)-1:nx*(i+1))
        ];
    end

    %% for scenario 4
    for i = 1:1:N
        ceq = [
            ceq;
            A*Xpred4(nx*i-1:nx*i)+B*Vvec4(i) + d4(:,t+i) - Xpred4(nx*(i+1)-1:nx*(i+1))
        ];
    end

    %% for scenario 5
    for i = 1:1:N
        ceq = [
            ceq;
            A*Xpred5(nx*i-1:nx*i)+B*Vvec5(i) + d5(:,t+i) - Xpred5(nx*(i+1)-1:nx*(i+1))
        ];
    end

    %% output 
    c = [];
end