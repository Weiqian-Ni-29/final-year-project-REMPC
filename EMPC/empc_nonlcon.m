function [c,ceq] = empc_nonlcon(s)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    global t N V A B;
    
    nx = 2;
    
    Xesti = s(1:nx*(t+1));
    Xpred = s(nx*t+1:nx*(t+N+1));
    
    Vvec = s(nx*(1+t+N)+1:nx*(t+N+1)+N);
    
    Gammavec = s(nx*(t+N+1)+N+1:end);
    
    % prediction
    ceq = [];
    for i = 1:1:N
        ceq = [ceq; A*Xpred(nx*i-1:nx*i)+B*Vvec(i) - Xpred(nx*(i+1)-1:nx*(i+1))];
    end
    
    % estimation
    if t>0
        for i = 1:1:t
            ceq = [ceq; A*Xesti(nx*i-1:nx*i)+B*V(i) - Xesti(nx*(i+1)-1:nx*(i+1))];
        end
    end
    
    
    % scaling factor
    c = [];
    
    for i = 1:N
        c = [c; min(1, Gammavec(i)+0.3) - Gammavec(i+1)];
    end
    
%     Rsettemp = A*Rset+B*(Vvec(1)+K*Xi)+Bw*Pw;
%     
%     for i = 1:1:N
%         Rsettempvertex = Rsettemp.V';
%         for j = 1:1:size(Rsettempvertex,2)
%             c = [c; Omega.H(:,1:nx)*(Rsettempvertex(:,j)-Xpred(nx*i+1:nx*i+nx))-Gammavec(i+1)*Omega.H(:,nx+1)];
%         end
%         Rsettemp = A*Rsettemp+B*(Vvec(i+1)+K*Xi)+Bw*Pw;
%     end
  
end