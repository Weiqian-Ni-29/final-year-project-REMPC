function [c,ceq] = empc_nonlcon(s)
    global t N V A B W_estimate begin_signal;
    
    nx = 2;
    Xpred = s(1:nx*(N+1)); % states for future timesteps (N+1)
    Vvec = s(nx*(1+N)+1:nx*(N+1)+N); % inputs for future timesteps (N)
   % load('multidisturbances.mat')
    W_est = zeros(2,N);
    if(t > 1)
        W_est(1,:) = W_estimate(1,t);
        W_est(2,:) = W_estimate(2,t);
    end
    W_est = W_est * 0.1 * begin_signal;
    % prediction
    ceq = [];
    for i = 1:1:N
        ceq = [ 
            ceq;
            A*Xpred(nx*i-1:nx*i)+B*Vvec(i) + W_est(:,i) - Xpred(nx*(i+1)-1:nx*(i+1)) % the difference in states between future timesteps and prediction
            ];
    end
    c = [];
end