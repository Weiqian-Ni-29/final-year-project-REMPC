function [fval] = empc_obj(s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    global t N Omegavertex coeff;
    
    nx = 2;
    
    Xpred = s(nx*t+1:nx*(t+N+1));
    
    Gammavec = s(nx*(t+N+1)+N+1:end);
    
%     fval = 0; % storage function
% 
%     for i = 0:N-1
%         Xtemp = Xpred(nx*i+1:nx*i+nx) + Gammavec(i+1)*Omegavertex;
%         fvaltemp = zeros(1,size(Omegavertex,2));
%         for j = 1:size(Omegavertex,2)
%             if Xtemp(2,j)<-2
%                 fvaltemp(j) = Xtemp(1,j) + 10*(Xtemp(2,j)+2)^2;
%             elseif Xtemp(2,j)>2
%                 fvaltemp(j) = Xtemp(1,j) + 10*(Xtemp(2,j)-2)^2;
%             else
%                 fvaltemp(j) = Xtemp(1,j) + 0;
%             end
%         end
%         
%         fval = fval + max(fvaltemp);
%     end
%     
    

%     
%     fval = -10000;
%     
%     for j = 1:size(Omegavertex,2)
%         fval = max(fval,Xpred(1)+Gammavec(1)*Omegavertex(1,j) + 0.01*(Xpred(nx)+Gammavec(1)*Omegavertex(2,j))^2);
%     end
    fval = 0;
    for i = 0:N-1
        f = 50*Xpred(nx*i+1);
        for j = 1:size(Omegavertex,2)
            f = max(f,50*(Xpred(nx*i+1)+Gammavec(i+1)*Omegavertex(1,j)));
        end
        if f == 50*Xpred(nx*i+1)
            error('problem')
        end
        fval = fval + f;
    end
    
    fval = fval + coeff(1)*Xpred(1)^2+coeff(2)*Xpred(1).*Xpred(2)+coeff(3)*Xpred(2).^2+coeff(4)*Xpred(1)+coeff(5)*Xpred(2);
%     fval = fval + (Xpred(nx*N+1:nx*N+nx))'*Q*(Xpred(nx*N+1:nx*N+nx));
end