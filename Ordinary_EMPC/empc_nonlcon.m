function [c,ceq] = empc_nonlcon(s)
    global t N V A B w1 w2;
    
    nx = 2;
    Xpred = s(1:nx*(N+1)); % states for future timesteps (N+1)
    Vvec = s(nx*(1+N)+1:nx*(N+1)+N); % inputs for future timesteps (N)
    load('multidisturbances.mat')
    % prediction
    ceq = [];
    for i = 1:1:N
        ceq = [ 
            ceq;
            A*Xpred(nx*i-1:nx*i)+B*Vvec(i) + w1*d1(:,t+i) + w2*d2(:,t+1) - Xpred(nx*(i+1)-1:nx*(i+1)) % the difference in states between future timesteps and prediction
            ];
    end
    c = [];
end