
function [ConvexProblem] = ParseSCvx(SCP,System,Reference,Cost,DiscreteDynamics,DiscreteConstraints,ScalingMatrices,dims,MosekCodes)

[Adynamics,bdynamics] = ParseDynamics(SCP,DiscreteDynamics,ScalingMatrices,dims);
dims.dyn = size(Adynamics,2);
[Aboundary,bboundary] = ParseBoundaryConditions(DiscreteConstraints,ScalingMatrices,dims);
[Alcs,blcs] = ParseLinearConstraints(DiscreteConstraints,ScalingMatrices,dims);
[At,btlc,btuc] = ParseTimeConstraints(SCP,System,ScalingMatrices,dims);

[Fsocs,gsocs,accsocs] = ParseSOCConstraints(DiscreteConstraints,ScalingMatrices,dims,MosekCodes);
[Fexps,gexps,accexps] = ParseExpConstraints(DiscreteConstraints,ScalingMatrices,dims,MosekCodes);
[Atr,btr,Ftr,gtr,acctr] = ParseTrustRegion(SCP,Reference,ScalingMatrices,dims,MosekCodes);
dims.n1 = size(Atr,2);
[c,cf,Fcost,gcost,acccost] = ParseCost(SCP,Cost,ScalingMatrices,dims,MosekCodes);

NN = length(c); % Number of Variables

% Linear Constraints
ConvexProblem.a = sparse([Adynamics zeros(size(Adynamics,1),NN-size(Adynamics,2));
                          Aboundary zeros(size(Aboundary,1),NN-size(Aboundary,2));
                          Alcs zeros(size(Alcs,1),NN-size(Alcs,2));
                          At zeros(size(At,1),NN-size(At,2));
                          Atr zeros(size(Atr,1),NN-size(Atr,2))]);

ConvexProblem.blc = [bdynamics;
                     bboundary;
                     -inf*ones(length(blcs),1);
                     btlc;
                     -inf*ones(length(btr),1)];

ConvexProblem.buc = [bdynamics;
                     bboundary;
                     blcs;
                     btuc;
                     btr];

% Cone Constraints:
ConvexProblem.f = sparse([Fsocs zeros(size(Fsocs,1),NN-size(Fsocs,2));
                          Fexps zeros(size(Fexps,1),NN-size(Fexps,2));
                          Ftr zeros(size(Ftr,1),NN-size(Ftr,2));
                          Fcost zeros(size(Fcost,1),NN-size(Fcost,2))]);

ConvexProblem.g = [gsocs;
                   gexps;
                   gtr;
                   gcost];

ConvexProblem.accs = [accsocs,...
                      accexps,...
                      acctr,...
                      acccost];
% ConvexProblem.f = sparse([Fsocs zeros(size(Fsocs,1),NN-size(Fsocs,2));
%                           Fcost zeros(size(Fcost,1),NN-size(Fcost,2))]);
% 
% ConvexProblem.g = [gsocs;
%                    gcost];
% 
% ConvexProblem.accs = [accsocs,...
%                       acccost];

% Cost:
ConvexProblem.c = c;
% ConvexProblem.c = zeros(length(c),1);
ConvexProblem.cf = cf;

end

function [A,b] = ParseDynamics(SCP,DiscreteDynamics,ScalingMatrices,dims)

% A Matrix:
A = zeros(dims.nx*(dims.N-1),dims.nx*dims.N);
for i = 1:dims.N-1
    A((i-1)*dims.nx+1:i*dims.nx,(i-1)*dims.nx+1:i*dims.nx) = DiscreteDynamics.A(:,:,i)*ScalingMatrices.Sx;
    A((i-1)*dims.nx+1:i*dims.nx,(i)*dims.nx+1:(i+1)*dims.nx) = -ScalingMatrices.Sx;
end

% B Matrix:
B = zeros(dims.nx*(dims.N-1),dims.nu*dims.N);
for i = 1:dims.N-1
    B((i-1)*dims.nx+1:i*dims.nx,(i-1)*dims.nu+1:i*dims.nu) = DiscreteDynamics.Bm(:,:,i)*ScalingMatrices.Su;
    B((i-1)*dims.nx+1:i*dims.nx,(i)*dims.nu+1:(i+1)*dims.nu) = DiscreteDynamics.Bp(:,:,i)*ScalingMatrices.Su;
end
A = [A B];

if SCP.AdaptiveMesh
    % H Matrix:
    H = zeros(dims.nx*(dims.N-1),dims.N);
    for i = 1:dims.N-1
        H((i-1)*dims.nx+1:i*dims.nx,i) = DiscreteDynamics.H(:,:,i)*ScalingMatrices.Sh;
    end
    A = [A H];

    % R vector:
    R = zeros((dims.N-1)*dims.nx,1);
    for i = 1:dims.N-1
        R((i-1)*dims.nx+1:i*dims.nx,1) = ScalingMatrices.cx-DiscreteDynamics.r(:,i)-DiscreteDynamics.A(:,:,i)*ScalingMatrices.cx-DiscreteDynamics.B(:,:,i)*ScalingMatrices.cu-DiscreteDynamics.H(:,:,i)*ScalingMatrices.ch;
    end
elseif SCP.FreeTime && ~SCP.AdaptiveMesh
    % F Matrix:
    F = zeros(dims.nx*(dims.N-1),1);
    for i = 1:dims.N-1
        F((i-1)*dims.nx+1:i*dims.nx) = DiscreteDynamics.F(:,i)*ScalingMatrices.Sp;
    end
    A = [A F];

    % R vector:
    R = zeros((dims.N-1)*dims.nx,1);
    for i = 1:dims.N-1
        R((i-1)*dims.nx+1:i*dims.nx,1) = ScalingMatrices.cx-DiscreteDynamics.r(:,i)-DiscreteDynamics.A(:,:,i)*ScalingMatrices.cx-DiscreteDynamics.Bm(:,:,i)*ScalingMatrices.cu-DiscreteDynamics.Bp(:,:,i)*ScalingMatrices.cu-DiscreteDynamics.F(:,i)*ScalingMatrices.cp;
    end
else

    R = zeros((dims.N-1)*dims.nx,1);
    for i = 1:dims.N-1
        R((i-1)*dims.nx+1:i*dims.nx,1) = ScalingMatrices.cx-DiscreteDynamics.r(:,i)-DiscreteDynamics.A(:,:,i)*ScalingMatrices.cx-DiscreteDynamics.Bm(:,:,i)*ScalingMatrices.cu-DiscreteDynamics.Bp(:,:,i)*ScalingMatrices.cu;
    end

end

% E Matrix: 
E = kron(eye(dims.N-1,dims.N-1),DiscreteDynamics.E);
A = [A E];

b = R;

end

function [A,b] = ParseBoundaryConditions(DiscreteConstraints,ScalingMatrices,dims)

A = [DiscreteConstraints.H0*ScalingMatrices.Sx zeros(dims.nic,dims.nx*(dims.N-1)) zeros(dims.nic,dims.dyn-dims.nx*dims.N+dims.ns) eye(dims.nic) zeros(dims.nic,dims.ntc)];
A = [A;
     zeros(dims.ntc,dims.nx*(dims.N-1)) DiscreteConstraints.Hf*ScalingMatrices.Sx zeros(dims.ntc,dims.dyn-dims.nx*dims.N+dims.ns) zeros(dims.ntc,dims.nic) eye(dims.ntc)];

b = [-DiscreteConstraints.l0 - DiscreteConstraints.H0*ScalingMatrices.cx;
     -DiscreteConstraints.lf - DiscreteConstraints.Hf*ScalingMatrices.cx];

end

function [A,b] = ParseLinearConstraints(DiscreteConstraints,ScalingMatrices,dims)

if isfield(DiscreteConstraints,'Kx')
    nl = size(DiscreteConstraints.Kx(:,:,1),1);
    Kx = zeros(nl*dims.N,dims.nx*dims.N);
    Ku = zeros(nl*dims.N,dims.nu*dims.N);
    lv = zeros(nl*dims.N,1);
    for i = 1:dims.N
        Kx((i-1)*nl+1:i*nl,(i-1)*dims.nx+1:i*dims.nx) = DiscreteConstraints.Kx(:,:,i)*ScalingMatrices.Sx;
        Ku((i-1)*nl+1:i*nl,(i-1)*dims.nu+1:i*dims.nu) = DiscreteConstraints.Ku(:,:,i)*ScalingMatrices.Su;
        lv((i-1)*nl+1:i*nl) = -DiscreteConstraints.l(:,i)-DiscreteConstraints.Kx(:,:,i)*ScalingMatrices.cx-DiscreteConstraints.Ku(:,:,i)*ScalingMatrices.cu;
    end
    
    A = [-Kx -Ku zeros(nl*dims.N,dims.dyn-dims.nx*dims.N-dims.nu*dims.N)];
    b = -lv;
else
    nl = 0;
    A = zeros(0,dims.dyn);
    b = [];
end

if dims.ns ~= 0

    C = zeros(dims.ns*dims.N,dims.nx*dims.N);
    D = zeros(dims.ns*dims.N,dims.nu*dims.N);
    rpv = zeros(dims.ns*dims.N);
    for i = 1:dims.N
        C((i-1)*dims.ns+1:i*dims.ns,(i-1)*dims.nx+1:i*dims.nx) = DiscreteConstraints.C(:,:,i)*ScalingMatrices.Sx;
        D((i-1)*dims.ns+1:i*dims.ns,(i-1)*dims.nu+1:i*dims.nu) = DiscreteConstraints.D(:,:,i)*ScalingMatrices.Su;
        rpv((i-1)*dims.ns+1:i*dims.ns) = -DiscreteConstraints.rp(:,i) - DiscreteConstraints.C(:,:,i)*ScalingMatrices.cx - DiscreteConstraints.D(:,:,i)*ScalingMatrices.cu;
    end

    A = [A zeros(nl*dims.N,dims.ns*dims.N);
         C D zeros(nl*dims.N,dims.dyn-dims.nx*dims.N-dims.nu*dims.N) eye(dims.ns*dims.N)];

    b = [b;
         rpv];

end

end

function [A,blc,buc] = ParseTimeConstraints(SCP,System,ScalingMatrices,dims)

if SCP.FreeTime && ~SCP.AdaptiveMesh
    A = [zeros(1,dims.nx*dims.N) zeros(1,dims.nu*dims.N) ScalingMatrices.Sp];
    blc = System.tmin - ScalingMatrices.cp;
    buc = System.tmax - ScalingMatrices.cp;
elseif SCP.FreeTime && SCP.AdaptiveMesh
    A = [zeros(1,dims.nx*dims.N) zeros(1,dims.nu*dims.N) ones(1,dims.N)*ScalingMatrices.Sh];
    blc = System.tmin - (dims.N-1)*ScalingMatrices.ch;
    buc = System.tmax - (dims.N-1)*ScalingMatrices.ch;
elseif ~SCP.FreeTime && SCP.AdaptiveMesh
    A = [zeros(1,dims.nx*dims.N) zeros(1,dims.nu*dims.N) ones(1,dims.N)*ScalingMatrices.Sh];
    blc = System.T - (dims.N-1)*ScalingMatrices.ch;
    buc = System.T - (dims.N-1)*ScalingMatrices.ch;
else
    A = [];
    blc = [];
    buc = [];
end

end

function [F,g,accs] = ParseSOCConstraints(DiscreteConstraints,ScalingMatrices,dims,MosekCodes)

F = [];
g = [];
accs = [];
if ~isfield(DiscreteConstraints,'Qx')
    return;
end

for k = 1:length(DiscreteConstraints.Qx)
    nl = size(DiscreteConstraints.Qx{k}(:,:,1),1);
    Qx = zeros(nl*dims.N,dims.nx*dims.N);
    Qu = zeros(nl*dims.N,dims.nu*dims.N);
    rv = zeros(nl*dims.N,1);
    for i = 1:dims.N
        Qx((i-1)*nl+1:i*nl,(i-1)*dims.nx+1:i*dims.nx) = DiscreteConstraints.Qx{k}(:,:,i)*ScalingMatrices.Sx;
        Qu((i-1)*nl+1:i*nl,(i-1)*dims.nu+1:i*dims.nu) = DiscreteConstraints.Qu{k}(:,:,i)*ScalingMatrices.Su;
        rv((i-1)*nl+1:i*nl) = DiscreteConstraints.r{k}(:,i)+DiscreteConstraints.Qx{k}(:,:,i)*ScalingMatrices.cx+DiscreteConstraints.Qu{k}(:,:,i)*ScalingMatrices.cu;
    end
    F = [F;
         Qx Qu zeros(nl*dims.N,dims.dyn-dims.nx*dims.N-dims.nu*dims.N)];
    g = [g;
         rv];
    accs = [accs kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE nl])];
end

end

function [F,g,accs] = ParseExpConstraints(DiscreteConstraints,ScalingMatrices,dims,MosekCodes)

F = [];
g = [];
accs = [];
if ~isfield(DiscreteConstraints,'Ex')
    return;
end

for k = 1:length(DiscreteConstraints.Ex)
    nl = size(DiscreteConstraints.Ex{k}(:,:,1),1);
    Ex = zeros(nl*dims.N,dims.nx*dims.N);
    Eu = zeros(nl*dims.N,dims.nu*dims.N);
    fv = zeros(nl*dims.N,1);
    for i = 1:dims.N
        Ex((i-1)*nl+1:i*nl,(i-1)*dims.nx+1:i*dims.nx) = DiscreteConstraints.Ex{k}(:,:,i)*ScalingMatrices.Sx;
        Eu((i-1)*nl+1:i*nl,(i-1)*dims.nu+1:i*dims.nu) = DiscreteConstraints.Eu{k}(:,:,i)*ScalingMatrices.Su;
        fv((i-1)*nl+1:i*nl) = DiscreteConstraints.f{k}(:,i)+DiscreteConstraints.Ex{k}(:,:,i)*ScalingMatrices.cx+DiscreteConstraints.Eu{k}(:,:,i)*ScalingMatrices.cu;
    end
    F = [F;
         Ex Eu zeros(nl*dims.N,dims.dyn-dims.nx*dims.N-dims.nu*dims.N)];
    g = [g;
         fv];
    accs = [accs kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_PRIMAL_EXP_CONE nl])];
end

end

function [Atr,btr,Ftr,gtr,accs] = ParseTrustRegion(SCP,Reference,ScalingMatrices,dims,MosekCodes)

if SCP.FreeTime && ~SCP.AdaptiveMesh
    Atr = [zeros(dims.N,dims.dyn+dims.ns+dims.nic+dims.ntc) SCP.alphax*eye(dims.N) SCP.alphau*eye(dims.N) SCP.alphap*ones(dims.N,1)];
    btr = ones(dims.N,1)*SCP.eta;
    
    Ftr = [kron(eye(dims.N),[zeros(1,dims.nx);eye(dims.nx)]) zeros(dims.N*(dims.nx+1),dims.dyn-dims.nx*dims.N+dims.ns+dims.nic+dims.ntc) kron(eye(dims.N),[1;zeros(dims.nx,1)]) zeros(dims.N*(dims.nx+1),dims.N) zeros(dims.N*(dims.nx+1),1);
           zeros(dims.N*(dims.nu+1),dims.N*dims.nx) kron(eye(dims.N),[zeros(1,dims.nu);eye(dims.nu)]) zeros(dims.N*(dims.nu+1),dims.dyn-dims.nx*dims.N-dims.nu*dims.N+dims.ns+dims.nic+dims.ntc) zeros(dims.N*(dims.nu+1),dims.N) kron(eye(dims.N),[1;zeros(dims.nu,1)]) zeros(dims.N*(dims.nu+1),1);
           zeros(2,dims.N*dims.nx) zeros(2,dims.N*dims.nu) [0;1] zeros(2,dims.dyn-dims.nx*dims.N-dims.nu*dims.N-1+dims.ns+dims.nic+dims.ntc) zeros(2,dims.N) zeros(2,dims.N) [1;0]];
    gtr = [-reshape([zeros(1,dims.N);ScalingMatrices.Sx^-1*(Reference.x-ScalingMatrices.cx)],[(dims.nx+1)*dims.N,1]);
           -reshape([zeros(1,dims.N);ScalingMatrices.Su^-1*(Reference.u-ScalingMatrices.cu)],[(dims.nu+1)*dims.N,1]);
           -[0;(ScalingMatrices.Sp^-1*(Reference.p-ScalingMatrices.cp))]];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nx+1]) kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nu+1]) MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 2];
elseif SCP.AdaptiveMesh
    Atr = [zeros(dims.N,dims.dyn+dims.ns+dims.nic+dims.ntc) SCP.alphax*eye(dims.N) SCP.alphau*eye(dims.N) SCP.alphah*eye(dims.N)];
    btr = ones(dims.N,1)*SCP.eta;

    Ftr = [kron(eye(dims.N),[zeros(1,dims.nx);eye(dims.nx)]) zeros(dims.N*(dims.nx+1),dims.dyn-dims.nx*dims.N+dims.ns+dims.nic+dims.ntc) kron(eye(dims.N),[1;zeros(dims.nx,1)]) zeros(dims.N*(dims.nx+1),dims.N) zeros(dims.N*(dims.nx+1),dims.N);
           zeros(dims.N*(dims.nu+1),dims.N*dims.nx) kron(eye(dims.N),[zeros(1,dims.nu);eye(dims.nu)]) zeros(dims.N*(dims.nu+1),dims.dyn-dims.nx*dims.N-dims.nu*dims.N+dims.ns+dims.nic+dims.ntc) zeros(dims.N*(dims.nu+1),dims.nx*dims.N) kron(eye(dims.N),[1;zeros(dims.nu,1)]) zeros(dims.N*(dims.nx+1),dims.N);
           zeros(2*(dims.N-1),dims.N*dims.nx) zeros(2*(dims.N-1),dims.N*dims.nu) kron(eye(dims.N-1),[0;1]) zeros(2*(dims.N-1),dims.dyn-dims.nx*dims.N-dims.nu*dims.N-dims.N+dims.ns+dims.nic+dims.ntc) zeros(2*(dims.N-1),dims.nx*dims.N) zeros(2*(dims.N-1),dims.nu*dims.N) kron(eye(dims.N-1),[1;0])];
    gtr = [-reshape([zeros(1,dims.N);ScalingMatrices.Sx^-1*(Reference.x-ScalingMatrices.cx)],[(dims.nx+1)*dims.N,1]);
           -reshape([zeros(1,dims.N);ScalingMatrices.Su^-1*(Reference.u-ScalingMatrices.cu)],[(dims.nu+1)*dims.N,1]);
           -reshape(([zeros(1,dims.N-1); ScalingMatrices.Sh^-1*(Reference.h-ScalingMatrices.ch)]),[2*(dims.N-1),1])];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nx+1]) kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nu+1]) kron(ones(1,dims.N-1),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 2])];
else
    Atr = [zeros(dims.N,dims.dyn+dims.ns+dims.nic+dims.ntc) SCP.alphax*eye(dims.N) SCP.alphau*eye(dims.N)];
    btr = ones(dims.N,1)*SCP.eta;

    Ftr = [kron(eye(dims.N),[zeros(1,dims.nx);eye(dims.nx)]) zeros(dims.N*(dims.nx+1),dims.dyn-dims.nx*dims.N+dims.ns+dims.nic+dims.ntc) kron(eye(dims.N),[1;zeros(dims.nx,1)]) zeros(dims.N*(dims.nx+1),dims.N);
           zeros(dims.N*(dims.nu+1),dims.N*dims.nx) kron(eye(dims.N),[zeros(1,dims.nu);eye(dims.nu)]) zeros(dims.N*(dims.nu+1),dims.dyn-dims.nx*dims.N-dims.nu*dims.N+dims.ns+dims.nic+dims.ntc) zeros(dims.N*(dims.nu+1),dims.N) kron(eye(dims.N),[1;zeros(dims.nu,1)])];
    gtr = [-reshape([zeros(1,dims.N);ScalingMatrices.Sx^-1*(Reference.x-ScalingMatrices.cx)],[(dims.nx+1)*dims.N,1]);
           -reshape([zeros(1,dims.N);ScalingMatrices.Su^-1*(Reference.u-ScalingMatrices.cu)],[(dims.nu+1)*dims.N,1])];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nx+1]) kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nu+1])];
end

end

function [c,cf,F,g,accs] = ParseCost(SCP,Cost,ScalingMatrices,dims,MosekCodes)

[lx,lu,ltf,Bx,Bu] = Cost();

c = [[zeros(dims.nx*(dims.N-1),1); ScalingMatrices.Sx*lx];
     [zeros(dims.nu*(dims.N-1),1); ScalingMatrices.Su*lu]];

cf = lx'*ScalingMatrices.cx + lu'*ScalingMatrices.cu;

if SCP.FreeTime && ~SCP.AdaptiveMesh
    c = [c;
         ltf*ScalingMatrices.Sp];
    cf = cf + ltf*ScalingMatrices.cp;
elseif SCP.FreeTime && SCP.AdaptiveMesh
    c = [c;
         ones(dims.N-1,1)*ScalingMatrices.Sh];
    cf = cf + ltf*ScalingMatrices.ch*(dims.N-1);
elseif ~SCP.FreeTime && ~SCP.AdaptiveMesh
    c = [c;
         0];
else
    c = [c;
         zeros(dims.N,1)];
end

if ~all(Bx(:)==0) && ~all(Bu(:)==0)
    F = [kron(eye(dims.N),[zeros(1,dims.nx);Bx*ScalingMatrices.Sx]) zeros(dims.N*(dims.nx+1),dims.n1-dims.nx*dims.N) kron(eye(dims.N),[1;zeros(dims.nx,1)])  zeros(dims.N*(dims.nx+1),dims.N);
          zeros(dims.N*(dims.nu+1),dims.N*dims.nx) kron(eye(dims.N),[zeros(1,dims.nu);Bu*ScalingMatrices.Su]) zeros(dims.N*(dims.nu+1),dims.n1-dims.nx*dims.N-dims.nu*dims.N) zeros(dims.N*(dims.nu+1),dims.N) kron(eye(dims.N),[1;zeros(dims.nu,1)])];
    g = [kron(ones(dims.N,1),[0;Bx*ScalingMatrices.cx]);
         kron(ones(dims.N,1),[0;Bu*ScalingMatrices.cu])];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 1+dims.nx]) kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 1+dims.nu])];

    c = [c;
         zeros(dims.n1-length(c),1);
         1/dims.N*ones(dims.N,1);
         1/dims.N*ones(dims.N,1)];

elseif ~all(Bx(:)==0) && all(Bu(:)==0)
    F = [kron(eye(dims.N),[zeros(1,dims.nx);Bx*ScalingMatrices.Sx]) zeros(dims.N*(dims.nx+1),dims.n1-dims.nx*dims.N) kron(eye(dims.N),[1;zeros(dims.nx,1)])];
    g = [kron(ones(dims.N,1),[0;Bx*ScalingMatrices.cx])];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 1+dims.nx])];

    c = [c;
         zeros(dims.n1-length(c),1);
         1/dims.N*ones(dims.N,1);];
elseif all(Bx(:)==0) && ~all(Bu(:)==0)
    F = [zeros(dims.N*(dims.nu+1),dims.N*dims.nx) kron(eye(dims.N),[zeros(1,dims.nu);Bu*ScalingMatrices.Su]) zeros(dims.N*(dims.nu+1),dims.n1-dims.nx*dims.N-dims.nu*dims.N) kron(eye(dims.N),[1;zeros(dims.nu,1)])];
    g = [kron(ones(dims.N,1),[0;Bu*ScalingMatrices.cu])];
    accs = [kron(ones(1,dims.N),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE 1+dims.nu])];

    c = [c;
         zeros(dims.n1-length(c),1);
         1/dims.N*ones(dims.N,1)];
else
    F = zeros(0,dims.n1);
    g = [];
    accs =[];

    c = [c;
         zeros(dims.n1-length(c),1)];
end

F = [F zeros(size(F,1),dims.N-1) zeros(size(F,1),logical(dims.ns)*dims.N) zeros(size(F,1),1) zeros(size(F,1),1);
     zeros((dims.nv+1)*(dims.N-1),dims.dyn-dims.nv*(dims.N-1)) kron(eye(dims.N-1),[zeros(1,dims.nv);eye(dims.nv)]) zeros((dims.nv+1)*(dims.N-1),size(F,2)-dims.dyn) kron(eye(dims.N-1),[1;zeros(dims.nv,1)]) zeros((dims.nv+1)*(dims.N-1),logical(dims.ns)*dims.N) zeros((dims.nv+1)*(dims.N-1),1) zeros((dims.nv+1)*(dims.N-1),1);
     zeros(logical(dims.ns)*(dims.ns+1)*dims.N,dims.dyn) kron(eye(dims.N),[zeros(1,dims.ns);eye(dims.ns)]) zeros(logical(dims.ns)*(dims.ns+1)*dims.N,size(F,2)-dims.dyn-dims.N*dims.ns) zeros(logical(dims.ns)*(dims.ns+1)*dims.N,dims.N-1) kron(eye(dims.N),[ones(logical(dims.ns),1);zeros(dims.ns,1)]) zeros(logical(dims.ns)*(dims.ns+1)*dims.N,1) zeros(logical(dims.ns)*(dims.ns+1)*dims.N,1);
     zeros(dims.nic+1,dims.dyn+dims.N*dims.ns) [zeros(1,dims.nic);eye(dims.nic)] zeros(dims.nic+1,size(F,2)-dims.dyn-dims.N*dims.ns-dims.nic) zeros(dims.nic+1,dims.N-1) zeros(dims.nic+1,logical(dims.ns)*dims.N) [1;zeros(dims.nic,1)] zeros(dims.nic+1,1);
     zeros(dims.ntc+1,dims.dyn+dims.N*dims.ns+dims.nic) [zeros(1,dims.ntc);eye(dims.ntc,dims.ntc)] zeros(dims.ntc+1,size(F,2)-dims.dyn-dims.N*dims.ns-dims.nic-dims.ntc) zeros(dims.ntc+1,dims.N-1) zeros(dims.ntc+1,logical(dims.ns)*dims.N) zeros(dims.ntc+1,1) [1;zeros(dims.ntc,1)]];
g = [g;
     zeros((dims.nv+1)*(dims.N-1),1);
     zeros(logical(dims.ns)*(dims.ns+1)*(dims.N),1);
     zeros(dims.nic+1,1);
     zeros(dims.ntc+1,1)];
accs = [accs kron(ones(1,dims.N-1),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nv+1]) kron(ones(1,dims.N*logical(dims.ns)),[MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.ns+1]) MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.nic+1 MosekCodes.symbcon.MSK_DOMAIN_QUADRATIC_CONE dims.ntc+1];

c = [c;
     SCP.lambda*ones(dims.N-1,1);
     SCP.lambdas*ones(dims.N*logical(dims.ns),1);
     SCP.lambdaic;
     SCP.lambdatc];

end