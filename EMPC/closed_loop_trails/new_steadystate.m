clear;
clc;
close all;

    global vertex gamma
    
    A = [1 1; 0 1];
    B = [0.5;1];
    xh = [8; 5];
    xl = [-8; -5];
    uh = 3;
    ul = -3;
    K = [-0.6696 -1.3261];
    
    Ax = [eye(2); -eye(2)];
    bx = [xh; -xl];
    Px = Polyhedron(Ax,bx);
    
    Pu = Polyhedron('V',[uh;ul]);
    
    load('Xi.mat')
    
    Pz = Px-Xi;
    Pv = Pu-K*Xi;
    
    zh = [7.5732; 4.725];
    zl = [-7.5732; -4.725];
    
    vh = 2.34956;
    vl = -2.34956;
    
    
    Aeq = [eye(2)-A -B];
    beq = zeros(2,1);
    
    lb = [zl;vl];
    ub = [zh;vh];
    
    s0 = zeros(3,1);
    
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunEvals',inf,'MaxIter',inf,'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-12);
%     s = fmincon(@(x)x(1),s0,[],[],Aeq,beq,lb,ub,[],options);
%     
%     zs = s(1:2)
%     vs = s(3)
%     save('SteadyState.mat','zs','vs')
    
    
    %%
    
%     load('Sets.mat','Omega');
    vertex = Xi.V';
    
    ZS = [];
    
    for gamma = [0.05:0.1:1 1]
        gamma
        
        [sopt, fval] = fmincon(@myfun,s0,[],[],Aeq,beq,lb,ub,[],options);

        zs = sopt(1:2)
        vs = sopt(3)
        fval

        f = zeros(1,size(vertex,2));

        for i = 1:1:size(vertex,2)
            x = [zs(1);zs(2)] + gamma*vertex(:,i);
            f(i) = x(1)+0.01*x(2)^2;
        end

        ind = find(f==fval);
        opti_ver = vertex(:,ind)
        ZS = [ZS zs];
    end
    
    
    z1 = zl(1):0.05:zh(1);
    z2 = zl(2):0.1:zh(2);
    [zz1,zz2] = meshgrid(z1,z2);
    zz = zz1 + 0.01*zz2.^2;
    surf(zz1,zz2,zz)
    
    save('newSteadyState.mat','zs','vs','opti_ver')
    
    
    function f = myfun(s)
        global vertex gamma
        
        
        f = s(1) + s(2)^2;
        
        for i = 1:1:size(vertex,2)
            x = [s(1);s(2)] + gamma*vertex(:,i);
            f = max(f, x(1)+0.01*x(2)^2);
        end
        
        if f == s(1) + 0.01*s(2)^2
            error('prob')
        end
        
    end
    
    