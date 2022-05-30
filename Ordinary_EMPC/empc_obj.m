function [fval] = empc_obj(s)
    global N;
    nx = 2;
    Xpred = s(nx+1:nx*(N+1));
    fval = 0;
    for i = 0:N-1
        f = Xpred(nx*i+1)^2 + (Xpred(nx*i+nx))^2;
        fval = fval + f;
    end
end