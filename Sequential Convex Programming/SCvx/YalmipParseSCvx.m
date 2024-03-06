
function [ConvexProblemOptimizer] = YalmipParseSCvx(SCP,System,ConvexConstraints,ScalingMatrices,dims)

yalmip('clear'); % Clear yalmip Variables

if ~SCP.RegenerateOptimizer && exist([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName '.mat'],'file') == 2
    fprintf('Saved Controller Detected, Loading...');
    load([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName]);
    fprintf('Done\n');
    return;
end

% Variables for Dynamic Matrices
A = sdpvar(dims.nx,dims.nx,dims.N-1,'full');
Bm = sdpvar(dims.nx,dims.nu,dims.N-1,'full');
Bp = sdpvar(dims.nx,dims.nu,dims.N-1,'full');
if SCP.FreeTime && ~SCP.AdaptiveMesh
    F = sdpvar(dims.nx,1,dims.N-1,'full');
end
if SCP.AdaptiveMesh
    H = sdpvar(dims.nx,1,dims.N-1,'full');
end
E = sdpvar(dims.nx,dims.nv,1,'full');
r = sdpvar(dims.nx,dims.N-1,'full');

% Variables for Cost
Aphi = sdpvar(1,dims.nx,'full');
Fphi = sdpvar(1,1,'full');
rphi = sdpvar(1,1,'full');
AL = sdpvar(1,dims.nx,dims.N,'full');
BL = sdpvar(1,dims.nu,dims.N,'full');
FL = sdpvar(1,1,dims.N,'full');
rL = sdpvar(1,dims.N,'full');

% Variables for Non Convex Constraint Matrices
if dims.ns ~= 0
    C = sdpvar(dims.ns,dims.nx,dims.N,'full');
    D = sdpvar(dims.ns,dims.nu,dims.N,'full');
    rp = sdpvar(dims.ns,dims.N,'full');
end

% Variables for Initial and Final Condition Constraint Matrices
H0 = sdpvar(dims.nic,dims.nx,'full');
Hf = sdpvar(dims.ntc,dims.nx,'full');
l0 = sdpvar(dims.nic,1,'full');
lf = sdpvar(dims.ntc,1,'full');

% Time Varying Problem Variables:
x = sdpvar(dims.nx,dims.N,'full'); % State
u = sdpvar(dims.nu,dims.N,'full'); % Input
if SCP.FreeTime && ~SCP.AdaptiveMesh
    p = sdpvar(1);
end
if SCP.AdaptiveMesh
    h = sdpvar(1,dims.N-1,'full');
end
c = sdpvar(dims.nc,1,'full');
v = sdpvar(dims.nv,dims.N-1,'full');
vs = sdpvar(dims.ns,dims.N,'full');
vic = sdpvar(dims.nic,1,'full');
vtc = sdpvar(dims.ntc,1,'full');

% Reference Solution Variables:
xbar = sdpvar(dims.nx,dims.N,'full'); % State
ubar = sdpvar(dims.nu,dims.N,'full'); % Input
if SCP.FreeTime && ~SCP.AdaptiveMesh
    pbar = sdpvar(1);
end
if SCP.AdaptiveMesh
    hbar = sdpvar(1,dims.N-1,'full');
end

% Trust Region Variable:
eta = sdpvar(1,1,'full');

% Constraints:
constraints = [];

% For Loop for Constraints
for k = 1:dims.N

    % Dynamics:
    if k < dims.N
        if ~SCP.FreeTime && ~SCP.AdaptiveMesh
            constraints = [constraints, ScalingMatrices.Sx*x(:,k+1)+ScalingMatrices.cx == A(:,:,k)*(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx) + Bm(:,:,k)*(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu) + Bp(:,:,k)*(ScalingMatrices.Su*u(:,k+1)+ScalingMatrices.cu) + E*v(:,k) + r(:,k)];
        end
        if SCP.FreeTime && ~SCP.AdaptiveMesh
             constraints = [constraints, ScalingMatrices.Sx*x(:,k+1)+ScalingMatrices.cx == A(:,:,k)*(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx) + Bm(:,:,k)*(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu) + Bp(:,:,k)*(ScalingMatrices.Su*u(:,k+1)+ScalingMatrices.cu) + F(:,1,k)*(ScalingMatrices.Sp*p+ScalingMatrices.cp) + E*v(:,k) + r(:,k)];
        end
        if SCP.AdaptiveMesh
            constraints = [constraints, ScalingMatrices.Sx*x(:,k+1)+ScalingMatrices.cx == A(:,:,k)*(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx) + Bm(:,:,k)*(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu) + Bp(:,:,k)*(ScalingMatrices.Su*u(:,k+1)+ScalingMatrices.cu) + H(:,1,k)*(ScalingMatrices.Sh*h(1,k)+ScalingMatrices.ch) + E*v(:,k) + r(:,k)];
        end
    end

    % NonConvex Constraints:
    if dims.ns ~= 0
        constraints = [constraints, vs(:,k) >= C(:,:,k)*(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx) + D(:,:,k)*(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu) + rp(:,k)];
    end

    % Convex Constraints: 
    constraints = [constraints, ConvexConstraints(System,ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx,ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu,c)];

    %Trust Region:
    if ~SCP.FreeTime && ~SCP.AdaptiveMesh
        constraints = [constraints, norm(x(:,k)-ScalingMatrices.Sx^-1*(xbar(:,k)-ScalingMatrices.cx),SCP.TrustRegionNorm) + norm(u(:,k)-ScalingMatrices.Su^-1*(ubar(:,k)-ScalingMatrices.cu),SCP.TrustRegionNorm)  <= eta];
    end
    if SCP.FreeTime && ~SCP.AdaptiveMesh
         % constraints = [constraints, norm(x(:,k)-ScalingMatrices.Sx^-1*(xbar(:,k)-ScalingMatrices.cx),SCP.TrustRegionNorm) + norm(u(:,k)-ScalingMatrices.Su^-1*(ubar(:,k)-ScalingMatrices.cu),SCP.TrustRegionNorm) + norm(p-ScalingMatrices.Sp^-1*(pbar-ScalingMatrices.cp),SCP.TrustRegionNorm) <= eta];
         constraints = [constraints, norm(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx-xbar(:,k),SCP.TrustRegionNorm) + norm(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu-ubar(:,k),SCP.TrustRegionNorm) + norm(ScalingMatrices.Sp*p+ScalingMatrices.cp-pbar,SCP.TrustRegionNorm)  <= eta];
    end
    if SCP.AdaptiveMesh
        constraints = [constraints, norm(x(:,k)-ScalingMatrices.Sx^-1*(xbar(:,k)-ScalingMatrices.cx),SCP.TrustRegionNorm) + norm(u(:,k)-ScalingMatrices.Su^-1*(ubar(:,k)-ScalingMatrices.cu),SCP.TrustRegionNorm) + norm(h(1,k)-ScalingMatrices.Sh^-1*(hbar(1,k)-ScalingMatrices.ch),SCP.TrustRegionNorm) <= eta];
    end

end

% Virtual Constraint:
constraints = [constraints, vs(:) >= 0];

% Boundary Conditions:
constraints = [constraints, 0 == H0*ScalingMatrices.Sx*x(:,1) + l0 + vic + H0*ScalingMatrices.cx];
constraints = [constraints, 0 == Hf*ScalingMatrices.Sx*x(:,dims.N) + lf + vtc + Hf*ScalingMatrices.cx];

% Time Constraints:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    T = ScalingMatrices.Sp*p + ScalingMatrices.cp;
    constraints = [constraints, System.tmin <= T <= System.tmax];
elseif SCP.FreeTime && SCP.AdaptiveMesh
    T = sum(ScalingMatrices.Sh*h+ScalingMatrices.ch);
    constraints = [constraints, System.tmin <= T <= System.tmax];
elseif ~SCP.FreeTime && SCP.AdaptiveMesh
    T = sum(ScalingMatrices.Sh*h+ScalingMatrices.ch);
    constraints = [constraints, System.T == T];
else
    T = System.T;
end

% Objective:
% Basic Terminal Cost:
phi = Aphi*(ScalingMatrices.Sx*x(:,dims.N)+ScalingMatrices.cx) + Fphi*(T) + rphi; 
% Augmented Terminal Cost:
phi = phi + SCP.lambdaic*norm(vic,1) + SCP.lambdatc*norm(vtc,1);

% Running Cost:
L = 0;
for k = 1:dims.N-1
    % Basic Running Cost:
    L = L + 1/dims.N*(AL(:,:,k)*(ScalingMatrices.Sx*x(:,k)+ScalingMatrices.cx) + BL(:,:,k)*(ScalingMatrices.Su*u(:,k)+ScalingMatrices.cu) + FL(:,:,k)*(T)+rL(:,k));
    % Augmented Running Cost:
    L = L + 1/dims.N*(SCP.lambda*norm(E*v(:,k),1));
    if dims.ns ~= 0
        L = L + 1/dims.N*(SCP.lambdas*sum(vs(:,k),1));
    end
end

% Full Objetive:
objective = phi + L;

% Parameters In and Solutions Out:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    parameters = {Aphi,Fphi,rphi,AL,BL,FL,rL,A,Bm,Bp,F,E,r,H0,Hf,l0,lf,xbar,ubar,pbar,eta};
    solutions = {objective,ScalingMatrices.Sx*x+ScalingMatrices.cx*ones(1,dims.N),ScalingMatrices.Su*u+ScalingMatrices.cu*ones(1,dims.N),ScalingMatrices.Sp*p+ScalingMatrices.cp,E*v,vic,vtc};
    if dims.nc ~= 0
        parameters = [parameters, {c}];
    end
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutions, {vs}];
    end
    
end
if SCP.FreeTime && SCP.AdaptiveMesh
    parameters = {Aphi,Fphi,rphi,AL,BL,FL,rL,A,Bm,Bp,H,E,r,H0,Hf,l0,lf,xbar,ubar,hbar,eta};
    solutions = {objective,ScalingMatrices.Sx*x+ScalingMatrices.cx*ones(1,dims.N),ScalingMatrices.Su*u+ScalingMatrices.cu*ones(1,dims.N),ScalingMatrices.Sh*h+ScalingMatrices.ch*ones(1,dims.N-1),E*v,vic,vtc};
    if dims.nc ~= 0
        parameters = [parameters, {c}];
    end
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutions, {vs}];
    end
end
if ~SCP.FreeTime
    parameters = {Aphi,Fphi,rphi,AL,BL,FL,rL,A,Bm,Bp,E,r,H0,Hf,l0,lf,xbar,ubar,eta};
    solutions = {objective,ScalingMatrices.Sx*x+ScalingMatrices.cx*ones(1,dims.N),ScalingMatrices.Su*u+ScalingMatrices.cu*ones(1,dims.N),E*v,vic,vtc};
    if dims.nc ~= 0
        parameters = [parameters, {c}];
    end
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutiobs, {vs}];
    end
end

ConvexProblemOptimizer = optimizer(constraints,objective,sdpsettings('solver','mosek'),parameters,solutions);

save([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName],'ConvexProblemOptimizer');

end