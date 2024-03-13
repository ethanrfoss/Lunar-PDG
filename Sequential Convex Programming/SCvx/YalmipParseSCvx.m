
function [ConvexProblemOptimizer,outFields] = YalmipParseSCvx(SCP,System,Convexities,ScalingMatrices,dims)

yalmip('clear'); % Clear yalmip Variables

% Check if optimizer needs to be regenerated
if ~SCP.RegenerateOptimizer && exist([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName '.mat'],'file') == 2
    fprintf('Saved Controller Detected, Loading...');
    load([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName]);
    fprintf('Done\n');
    return;
end

% Variables for Dynamic Matrices
Vars.A = NewVar([dims.nx,dims.nx,dims.N-1],'full');
Vars.Bm = NewVar([dims.nx,dims.nu,dims.N-1],'full');
Vars.Bp = NewVar([dims.nx,dims.nu,dims.N-1],'full');
Vars.C = NewVar([dims.nx,dims.np,dims.N-1],'full');
Vars.D = NewVar([dims.nx,dims.nh,dims.N-1],'full');
Vars.E = NewVar([dims.nx,dims.nv,1],'full');
Vars.r = NewVar([dims.nx,dims.N-1],'full');
Vars.G = NewVar([dims.nx,dims.nx,dims.N-1],'full');

% Variables for Cost
Vars.Aphi = NewVar([1,dims.nx],'full');
Vars.Bphi = NewVar([1,dims.nu],'full');
Vars.Cphi = NewVar([1,dims.np],'full');
Vars.Dphi = NewVar([1,dims.nh],'full');
Vars.rphi = NewVar([1,1],'full');
Vars.AL = NewVar([1,dims.nx,dims.N],'full');
Vars.BL = NewVar([1,dims.nu,dims.N],'full');
Vars.CL = NewVar([1,dims.np,dims.N],'full');
Vars.DL = NewVar([1,dims.nh,dims.N],'full');
Vars.rL = NewVar([1,dims.N],'full');

% Variables for Non Convex Constraint Matrices
Vars.As = NewVar([dims.ns,dims.nx,dims.N],'full');
Vars.Bs = NewVar([dims.ns,dims.nu,dims.N],'full');
Vars.Cs = NewVar([dims.ns,dims.np,dims.N],'full');
Vars.Ds = NewVar([dims.ns,dims.nh,dims.N],'full');
Vars.rs = NewVar([dims.ns,dims.N],'full');

% Variables for Boundary Condition Matrices
Vars.H0 = NewVar([dims.nic,dims.nx],'full');
Vars.I0 = NewVar([dims.nic,dims.nu],'full');
Vars.J0 = NewVar([dims.nic,dims.np],'full');
Vars.K0 = NewVar([dims.nic,dims.nh],'full');
Vars.l0 = NewVar([dims.nic,1],'full');

Vars.Hf = NewVar([dims.ntc,dims.nx],'full');
Vars.If = NewVar([dims.ntc,dims.nu],'full');
Vars.Jf = NewVar([dims.ntc,dims.np],'full');
Vars.Kf = NewVar([dims.ntc,dims.nh],'full');
Vars.lf = NewVar([dims.ntc,1],'full');

% Time Varying Problem Variables:
Vars.x = NewVar([dims.nx,dims.N],'full'); % State
Vars.u = NewVar([dims.nu,dims.N],'full'); % Input
Vars.p = NewVar([dims.np,1],'full'); % Parameter
Vars.h = NewVar([dims.nh,dims.N],'full'); % Time Varying Parameter
Vars.v = NewVar([dims.nv,dims.N-1],'full');
Vars.vs = NewVar([dims.ns,dims.N],'full');
Vars.vic = NewVar([dims.nic,1],'full');
Vars.vtc = NewVar([dims.ntc,1],'full');

% Reference Solution Variables:
Vars.xbar = NewVar([dims.nx,dims.N],'full'); % State
Vars.ubar = NewVar([dims.nu,dims.N],'full'); % Input
Vars.pbar = NewVar([dims.np,1],'full'); % Parameter
Vars.hbar = NewVar([dims.nh,dims.N],'full'); % Time Varying Parameter

% Continuation Parameter:
Vars.c = NewVar([dims.nc,1],'full');

% Trust Region Variable:
Vars.eta = NewVar([1,1],'full');

% Constraints:
constraints = [];

% For Loop for Constraints
for k = 1:dims.N

    % Dynamics:
    if k < dims.N
        constraints = [constraints, ScalingMatrices.Sx*Vars.x(:,k+1)+ScalingMatrices.cx == Vars.A(:,:,k)*(ScalingMatrices.Sx*Vars.x(:,k)+ScalingMatrices.cx) + Vars.Bm(:,:,k)*(ScalingMatrices.Su*Vars.u(:,k)+ScalingMatrices.cu) + Vars.Bp(:,:,k)*(ScalingMatrices.Su*Vars.u(:,k+1)+ScalingMatrices.cu) + Vars.C(:,:,k)*(ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp) + Vars.D(:,:,k)*(ScalingMatrices.Sh*Vars.h(:,k) + ScalingMatrices.ch) + Vars.E*Vars.v(:,k) + Vars.r(:,k)];
    end

    constraints = [constraints, Vars.vs(:,k) >= Vars.As(:,:,k)*(ScalingMatrices.Sx*Vars.x(:,k)+ScalingMatrices.cx) + Vars.Bs(:,:,k)*(ScalingMatrices.Su*Vars.u(:,k)+ScalingMatrices.cu) + Vars.Cs(:,:,k)*(ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp) + Vars.Ds(:,:,k)*(ScalingMatrices.Sh*Vars.h(:,k)+ScalingMatrices.ch) + Vars.rs(:,k)];

    %Trust Region:
    constraints = [constraints, SCP.alphax*norm(Vars.x(:,k)-ScalingMatrices.Sx^-1*(Vars.xbar(:,k)-ScalingMatrices.cx),SCP.TrustRegionNorm) + SCP.alphau*norm(Vars.u(:,k)-ScalingMatrices.Su^-1*(Vars.ubar(:,k)-ScalingMatrices.cu),SCP.TrustRegionNorm) + SCP.alphap*norm(Vars.p-ScalingMatrices.Sp^-1*(Vars.pbar-ScalingMatrices.cp),SCP.TrustRegionNorm) + SCP.alphah*norm(Vars.h(:,k)-ScalingMatrices.Sh^-1*(Vars.hbar(:,k)-ScalingMatrices.ch),SCP.TrustRegionNorm)  <= Vars.eta];

end

% Virtual Constraint:
constraints = [constraints, Vars.vs(:) >= 0];

% Boundary Conditions:
constraints = [constraints, 0 == Vars.H0*(ScalingMatrices.Sx*Vars.x(:,1) + ScalingMatrices.cx) + Vars.I0*(ScalingMatrices.Su*Vars.u(:,1) + ScalingMatrices.cu) + Vars.J0*(ScalingMatrices.Sp*Vars.p + ScalingMatrices.cp) + Vars.K0*(ScalingMatrices.Sh*Vars.h(:,1) + ScalingMatrices.ch) + Vars.l0 + Vars.vic];
constraints = [constraints, 0 == Vars.Hf*(ScalingMatrices.Sx*Vars.x(:,end) + ScalingMatrices.cx) + Vars.If*(ScalingMatrices.Su*Vars.u(:,end) + ScalingMatrices.cu) + Vars.Jf*(ScalingMatrices.Sp*Vars.p + ScalingMatrices.cp) + Vars.Kf*(ScalingMatrices.Sh*Vars.h(:,end) + ScalingMatrices.ch) + Vars.lf + Vars.vtc];

% Objective:
% Basic Terminal Cost:
phi = Vars.Aphi*(ScalingMatrices.Sx*Vars.x(:,dims.N)+ScalingMatrices.cx) + Vars.Bphi*(ScalingMatrices.Su*Vars.u(:,dims.N)+ScalingMatrices.cu) + Vars.Cphi*(ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp) + Vars.Dphi*(ScalingMatrices.Sh*Vars.h(:,dims.N)+ScalingMatrices.ch) + Vars.rphi; 
% Augmented Terminal Cost:
phi = phi + SCP.lambdaic*norm(Vars.vic,1) + SCP.lambdatc*norm(Vars.vtc,1);

% Running Cost:
L = 0;
for k = 1:dims.N
    % Basic Running Cost:
    L = L + 1/dims.N*(Vars.AL(:,:,k)*(ScalingMatrices.Sx*Vars.x(:,k)+ScalingMatrices.cx) + Vars.BL(:,:,k)*(ScalingMatrices.Su*Vars.u(:,k)+ScalingMatrices.cu) + Vars.CL(:,:,k)*(ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp) + Vars.DL(:,:,k)*(ScalingMatrices.Sh*Vars.h(:,k)+ScalingMatrices.ch) + Vars.rL(:,k));
    % Augmented Running Cost:
    if k < dims.N
        L = L + 1/dims.N*(SCP.lambda*norm(Vars.E*Vars.v(:,k),1));
    end
    L = L + 1/dims.N*(SCP.lambdas*sum(Vars.vs(:,k),1));
end

% Full Objetive:
objective = phi + L;

% Get User Specified Convexities:
[AdditionalConstraints,AdditionalObjective,Vars] = Convexities(System,ScalingMatrices.Sx*Vars.x+ScalingMatrices.cx*ones(1,dims.N),ScalingMatrices.Su*Vars.u+ScalingMatrices.cu*ones(1,dims.N),ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp,ScalingMatrices.Sh*Vars.h+ScalingMatrices.ch*ones(1,dims.N),Vars.c,Vars);
constraints = [constraints, AdditionalConstraints];
objective = objective + AdditionalObjective;

% Parameters In and Solutions Out:
parameters = {Vars.Aphi,Vars.Bphi,Vars.Cphi,Vars.Dphi,Vars.rphi,Vars.AL,Vars.BL,Vars.CL,Vars.DL,Vars.rL,Vars.A,Vars.Bm,Vars.Bp,Vars.C,Vars.D,Vars.E,Vars.r,Vars.G,Vars.H0,Vars.I0,Vars.J0,Vars.K0,Vars.l0,Vars.Hf,Vars.If,Vars.Jf,Vars.Kf,Vars.lf,Vars.As,Vars.Bs,Vars.Cs,Vars.Ds,Vars.rs,Vars.xbar,Vars.ubar,Vars.pbar,Vars.hbar,Vars.c,Vars.eta};
outVars = rmfield(Vars,{'x','u','p','h','c'});
solutions = [{objective,ScalingMatrices.Sx*Vars.x+ScalingMatrices.cx*ones(1,dims.N),ScalingMatrices.Su*Vars.u+ScalingMatrices.cu*ones(1,dims.N),ScalingMatrices.Sp*Vars.p+ScalingMatrices.cp,ScalingMatrices.Sh*Vars.h+ScalingMatrices.ch*ones(1,dims.N)},Vars.c,struct2cell(outVars)'];

% Remove Empty Parameters and Solutions:
parameters = parameters(~cellfun('isempty',parameters));
solutions = solutions(~cellfun('isempty',solutions));
VarNames = fieldnames(outVars);
outFields = VarNames(find(~cellfun('isempty',struct2cell(outVars))));

% Create Optimizer:
ConvexProblemOptimizer = optimizer(constraints,objective+AdditionalObjective,sdpsettings('solver','mosek'),parameters,solutions);

% Save Optimizer:
save([fileparts(mfilename('fullpath')) '\Optimizers\' SCP.OptimizerName],'ConvexProblemOptimizer');

end

function var = NewVar(dim,type)

if any(dim==0)
    var = zeros(dim);
else
    dimcell = num2cell(dim);
    var = sdpvar(dimcell{:},type);
end

end