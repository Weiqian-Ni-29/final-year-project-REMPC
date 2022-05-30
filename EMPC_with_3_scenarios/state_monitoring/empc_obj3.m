function [fval] = empc_obj3(s)
    % f = x1(t)^2+x2(t)^2
    global N w1 w2 w3;
    nx = 2;

    %% scenario 1
    Xpred1 = s(nx+1:nx*(N+1));  % predicted states for the first scenario
    fval1 = 0;
    for i = 0:N-1
        f = Xpred1(nx*i+1)^2 + (Xpred1(nx*i+nx))^2;
        fval1 = fval1 + f;
    end

    %% scenario 2
    Xpred2 = s(nx*(N+1)+N+nx+1:2*nx*(N+1)+N); % predicted states for the second scenario
    fval2 = 0;
    for i = 0:N-1
        f = Xpred2(nx*i+1)^2 + (Xpred2(nx*i+nx))^2;
        fval2 = fval2 + f;
    end

    %% scenario 3
    Xpred3 = s(2*(nx*(N+1)+N)+nx+1:2*(nx*(N+1)+N)+nx*(N+1)); % predicted states for the second scenario
    fval3 = 0;
    for i = 0:N-1
        f = Xpred3(nx*i+1)^2 + (Xpred3(nx*i+nx))^2;
        fval3 = fval3 + f;
    end

    %% compute the overall cost as a linear combination of both
    fval = w1*fval1+w2*fval2+w3*fval3;
end