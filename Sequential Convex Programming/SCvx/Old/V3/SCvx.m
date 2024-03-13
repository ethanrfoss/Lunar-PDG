
function [L,Solution] = SCvx(SCP,System,InitialGuess,Variables,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition,Scaling,Continuation)

% Algorithm Start
tic;

% Determine Problem Dimensions
[dims] = Dimensions(System,SCP,InitialGuess,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling);

% Initialize Continuation Parameter:
if dims.nc == 0
    gamma = 1;
else
    gamma = 0;
    InitialGuess.c = Continuation(System,gamma);
end

% Find Jacobians
fprintf('Linearizing Problem...');
[LinearizedDynamics,LinearizedConstraints] = Linearize(SCP,System,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition);
fprintf('Done\n');

% Parse the Problem Using YALMIP
fprintf('Parsing Problem...');
[Controller] = Parse(SCP,System,Cost,ConvexConstraints,Scaling,dims);
fprintf('Done\n');

% Discretize
fprintf('Discretizing Problem...');
[DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,InitialGuess,LinearizedDynamics,LinearizedConstraints,dims);
fprintf('Done\n');

% Nonlinear Cost:
Jbar = NonLinearCost(SCP,InitialGuess,LinearizedConstraints,Cost,delta,dims);

Reference = InitialGuess;
iter = 1;

while true
    
    % Solve Problem as a SOCP
    fprintf('Solving Convex Subproblem..');
    if SCP.FreeTime && ~SCP.AdaptiveMesh
        inputs = {DiscreteDynamics.A,...
                  DiscreteDynamics.Bm,...
                  DiscreteDynamics.Bp,...
                  DiscreteDynamics.F,...
                  DiscreteDynamics.E,...
                  DiscreteDynamics.r,...
                  DiscreteConstraints.H0,...
                  DiscreteConstraints.Hf,...
                  DiscreteConstraints.l0,...
                  DiscreteConstraints.lf,...
                  Reference.x,...
                  Reference.u,...
                  Reference.p,...
                  Reference.c,...
                  SCP.eta};
        if dims.ns ~= 0
            inputs = [inputs,...
                      {DiscreteConstraints.C,...
                      DiscreteConstraints.D,...
                      DiscreteConstraints.rp}];
        end
        [solution,infeasible] = Controller{inputs};
        if infeasible
            error('Convex Subproblem was infeasible!');
        end
        Lstar = solution{1};
        NewReference.x = solution{2};
        NewReference.u = solution{3};
        NewReference.p = solution{4};
        NewReference.v = solution{5};
        NewReference.vic = solution{6};
        NewReference.vtc = solution{7};
        NewReference.c = Reference.c;
    end
    if SCP.FreeTime && SCP.AdaptiveMesh
        inputs = {DiscreteDynamics.A,...
                  DiscreteDynamics.Bm,...
                  DiscreteDynamics.Bp,...
                  DiscreteDynamics.H,...
                  DiscreteDynamics.E,...
                  DiscreteDynamics.r,...
                  DiscreteConstraints.H0,...
                  DiscreteConstraints.Hf,...
                  DiscreteConstraints.l0,...
                  DiscreteConstraints.lf,...
                  Reference.x,...
                  Reference.u,...
                  Reference.h,...
                  Reference.c,...
                  SCP.eta};
        if dims.ns ~= 0
            inputs = [inputs,...
                      {DiscreteConstraints.C,...
                      DiscreteConstraints.D,...
                      DiscreteConstraints.rp}];
        end
        [solution,infeasible] = Controller{inputs};
        if infeasible
            error('Convex Subproblem was infeasible!');
        end
        Lstar = solution{1};
        NewReference.x = solution{2};
        NewReference.u = solution{3};
        NewReference.p = solution{4};
        NewReference.v = solution{5};
        NewReference.vic = solution{6};
        NewReference.vtc = solution{7};
        NewReference.c = Reference.c;
    end
    if ~SCP.FreeTime
        inputs = {DiscreteDynamics.A,...
                  DiscreteDynamics.Bm,...
                  DiscreteDynamics.Bp,...
                  DiscreteDynamics.E,...
                  DiscreteDynamics.r,...
                  DiscreteConstraints.H0,...
                  DiscreteConstraints.Hf,...
                  DiscreteConstraints.l0,...
                  DiscreteConstraints.lf,...
                  Reference.x,...
                  Reference.u,...
                  Reference.c,...
                  SCP.eta};
        if dims.ns ~= 0
            inputs = [inputs,...
                      {DiscreteConstraints.C,...
                      DiscreteConstraints.D,...
                      DiscreteConstraints.rp}];
        end
        [solution,infeasible] = Controller{inputs};
        if infeasible
            error('Convex Subproblem was infeasible!');
        end
        Lstar = solution{1};
        NewReference.x = solution{2};
        NewReference.u = solution{3};
        NewReference.p = solution{4};
        NewReference.v = solution{5};
        NewReference.vic = solution{6};
        NewReference.vtc = solution{7};
        NewReference.c = Reference.c;
    end
    fprintf('Done\n');

    % Check Convergence Criterion
    if abs(Jbar-Lstar) <= SCP.epsilon && gamma == 1 %norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)) <= SCP.epsilon
        disp(['SCP Algorithm Converged in ' num2str(toc) ' seconds']);
        break;
    else
        disp(['Current Iteration: ' num2str(iter) ' Current Tolerance: ' num2str(abs(Jbar-Lstar)) ' Current Trust Region: ' num2str(SCP.eta)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
    end

    if iter >= SCP.maxIter
        disp(['Maximum Iteration Count Reached. Current Tolerance: ' num2str(abs(Jbar-Lstar))]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
        break
    end
    
    % Redescretize with solution
    fprintf('Discretizing Problem...');
    [DiscreteDynamicsStar,DiscreteConstraintsStar,delta] = Discretize(SCP,NewReference,LinearizedDynamics,LinearizedConstraints,dims);
    fprintf('Done\n');

    % New Nonlinear Cost
    Jstar = NonLinearCost(SCP,NewReference,LinearizedConstraints,Cost,delta,dims);
    
    % Calculate Convexification Accuracy
    rho = (Jbar - Jstar)/(Jbar-Lstar);
    
    % Update trust region
     if rho < SCP.rho0
        % Update eta, dont update trajectory
        SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
    else
        if rho < SCP.rho1
            % Update eta
            SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
        end
        % If rho1 < rho <= rho2, do not update eta
        if rho >= SCP.rho2
            % Update eta
            SCP.eta = min(SCP.eta1,SCP.betagr*SCP.eta);
        end
        % Update Nominal Solution
        Reference = NewReference;
        % Update Nonlinear Cost
        Jbar = Jstar;
        % Update Discretiszation
        DiscreteDynamics = DiscreteDynamicsStar;
        DiscreteConstraints = DiscreteConstraintsStar;
    end
    % Update Result:
    L(iter) = Lstar;
    Solution(iter) = Reference;

    % Apply Continuation:
    if dims.nc ~= 0
        if max(vecnorm(delta)) <= SCP.gammaeps && gamma < 1
            gamma = min(gamma + SCP.dgamma,1);
            disp(['Applying Continuation. gamma = ' num2str(gamma)]);
            Reference.c = Continuation(System,gamma);
        end
    end

    iter = iter + 1;
    
end


end

function J = NonLinearCost(SCP,ReferenceSolution,LinearizedConstraints,Cost,delta,dims)

% Cost:
[lx,lu,ltf,Bx,Bu] = Cost();

% Terminal Cost:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    phi = lx'*ReferenceSolution.x(:,end) + lu'*ReferenceSolution.u(:,end) + ltf*ReferenceSolution.p + SCP.lambdaic*norm(LinearizedConstraints.g0(ReferenceSolution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(ReferenceSolution.x(:,end)),2);
end
if SCP.FreeTime && SCP.AdaptiveMesh
    phi = lx'*ReferenceSolution.x(:,end) + lu'*ReferenceSolution.u(:,end) + ltf*sum(ReferenceSolution.h) + SCP.lambdaic*norm(LinearizedConstraints.gic(ReferenceSolution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(ReferenceSolution.x(:,end)),2);
end
if ~SCP.FreeTime
    phi = lx'*ReferenceSolution.x(:,end) + lu'*ReferenceSolution.u(:,end) + SCP.lambdaic*norm(LinearizedDynamics.gic(ReferenceSolution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedDynamics.gf(ReferenceSolution.x(:,end)),2);
end

Gamma = 0;
% Running Cost:
for k = 1:dims.N-1
    if dims.ns ~= 0
        Gamma = Gamma + 1/dims.N*(norm(Bx*ReferenceSolution.x(:,k)) + norm(Bu*ReferenceSolution.u(:,k)) + SCP.lambda*norm(delta(:,k)) + SCP.lambdas*norm(max(LinearizedConstraints.s(ReferenceSolution.x(:,k),ReferenceSolution.u(:,k)),0)));
    else
        Gamma = Gamma + 1/dims.N*(norm(Bx*ReferenceSolution.x(:,k)) + norm(Bu*ReferenceSolution.u(:,k)) + SCP.lambda*norm(delta(:,k)));
    end
end
if dims.ns ~= 0
    Gamma = Gamma + 1/dims.N*(norm(Bx*ReferenceSolution.x(:,k)) + norm(Bu*ReferenceSolution.u(:,k)) + SCP.lambdas*norm(max(LinearizedConstraints.s(ReferenceSolution.x(:,k),ReferenceSolution.u(:,k)),0)));
else
    Gamma = Gamma + 1/dims.N*(norm(Bx*ReferenceSolution.x(:,k)) + norm(Bu*ReferenceSolution.u(:,k)));
end

% Total Cost:
J = phi + Gamma;

end

function [Controller] = Parse(SCP,System,Cost,ConvexConstraints,Scaling,dims)
tic;
yalmip('clear'); % Clear yalmip Variables

if ~SCP.RegenerateController && exist([fileparts(mfilename('fullpath')) '\Controllers\' SCP.ControllerName '.mat'],'file') == 2
    load([fileparts(mfilename('fullpath')) '\Controllers\' SCP.ControllerName '.mat'],'Controller');
    return;
end

% Get Scaling Matrices:
[Sx,Su,cx,cu] = ScalingMatrices(System,Scaling);
if SCP.FreeTime && ~SCP.AdaptiveMesh
    Sp = System.tmax-System.tmin;
    cp = System.tmin;
end
if SCP.AdaptiveMesh
    Sh = System.MaxStep-System.MinStep;
    ch = System.MinStep;
end

% Get Cost Parameters:
[lx,lu,ltf,Bx,Bu] = Cost();

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

% Variables for Non Convex Constraint Matrices
if dims.ns ~= 0
    C = sdpvar(dims.ns,dims.nx,dims.N,'full');
    D = sdpvar(dims.ns,dims.nu,dims.N,'full');
    rp = sdpvar(dims.ns,1,dims.N,'full');
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
            constraints = [constraints, Sx*x(:,k+1)+cx == A(:,:,k)*(Sx*x(:,k)+cx) + Bm(:,:,k)*(Su*u(:,k)+cu) + Bp(:,:,k)*(Su*u(:,k+1)+cu) + E*v(:,k) + r(:,k)];
        end
        if SCP.FreeTime && ~SCP.AdaptiveMesh
             constraints = [constraints, Sx*x(:,k+1)+cx == A(:,:,k)*(Sx*x(:,k)+cx) + Bm(:,:,k)*(Su*u(:,k)+cu) + Bp(:,:,k)*(Su*u(:,k+1)+cu) + F(:,1,k)*(Sp*p+cp) + E*v(:,k) + r(:,k)];
        end
        if SCP.AdaptiveMesh
            constraints = [constraints, Sx*x(:,k+1)+cx == A(:,:,k)*(Sx*x(:,k)+cx) + Bm(:,:,k)*(Su*u(:,k)+cu) + Bp(:,:,k)*(Su*u(:,k+1)+cu) + H(:,1,k)*(Sh*h(1,k)+ch) + E*v(:,k) + r(:,k)];
        end
    end

    % NonConvex Constraints:
    if dims.ns ~= 0
        constraints = [constraints, vs(:,k) >= C(:,:,k)*Sx*x(:,k) + D(:,:,k)*Su*u(:,k) + rp(:,k) + C(:,:,k)*cx + D(:,:,k)*cu];
    end

    % Convex Constraints: 
    constraints = [constraints, ConvexConstraints(System,Sx*x(:,k)+cx,Su*u(:,k)+cu,c)];

    %Trust Region:
    if ~SCP.FreeTime && ~SCP.AdaptiveMesh
        constraints = [constraints, norm(x(:,k)-Sx^-1*(xbar(:,k)-cx),SCP.TrustRegionNorm) + norm(u(:,k)-Su^-1*(ubar(:,k)-cu),SCP.TrustRegionNorm)  <= eta];
    end
    if SCP.FreeTime && ~SCP.AdaptiveMesh
         constraints = [constraints, norm(x(:,k)-Sx^-1*(xbar(:,k)-cx),SCP.TrustRegionNorm) + norm(u(:,k)-Su^-1*(ubar(:,k)-cu),SCP.TrustRegionNorm) + norm(p-Sp^-1*(pbar-cp),SCP.TrustRegionNorm) <= eta];
    end
    if SCP.AdaptiveMesh
        constraints = [constraints, norm(x(:,k)-Sx^-1*(xbar(:,k)-cx),SCP.TrustRegionNorm) + norm(u(:,k)-Su^-1*(ubar(:,k)-cu),SCP.TrustRegionNorm) + norm(h(1,k)-Sh^-1*(hbar(1,k)-ch),SCP.TrustRegionNorm) <= eta];
    end

end

% Boundary Conditions:
constraints = [constraints, 0 == H0*Sx*x(:,1) + l0 + vic + H0*cx];
constraints = [constraints, 0 == Hf*Sx*x(:,dims.N) + lf + vtc + Hf*cx];

% Time Constraints:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    constraints = [constraints, System.tmin <= Sp*p + cp <= System.tmax];
end
if SCP.FreeTime && SCP.AdaptiveMesh
    constraints = [constraints, System.tmin <= sum(Sh*h+ch) <= System.tmax];
end
if ~SCP.FreeTime && SCP.AdaptiveMesh
    constraints = [constraints, System.T == sum(Sh*h+ch)];
end

% Objective:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    objective = lx'*(Sx*x(:,dims.N)+cx) + lu'*(Su*u(:,dims.N)+cu) + ltf*(Sp*p + cp) + SCP.lambdaic*norm(vic,2) + SCP.lambdatc*norm(vtc,2);
end
if SCP.FreeTime && SCP.AdaptiveMesh
    objective = lx'*(Sx*x(:,dims.N)+cx) + lu'*(Su*u(:,dims.N)+cu) + ltf*sum(Sh*h+ch) + SCP.lambdaic*norm(vic,2) + SCP.lambdatc*norm(vtc,2);
end
if ~SCP.FreeTime
    objective = lx'*(Sx*x(:,dims.N)+cx) + lu'*(Su*u(:,dims.N)+cu) + SCP.lambdaic*norm(vic,2) + SCP.lambdatc*norm(vtc,2);
end

for k = 1:dims.N-1
    objective = objective + 1/dims.N*(norm(Bx*(Sx*x(:,k)+cx),2) + norm(Bu*(Su*u(:,k)+cu),2) + SCP.lambda*norm(E*v(:,k),2));
    if dims.ns ~= 0
        objective = objective + 1/dims.N(SCP.lambdas*norm(vs(:,k),2));
    end
end
objective = objective + 1/dims.N*(norm(Bx*(Sx*x(:,k)+cx),2)+norm(Bu*(Su*u(:,k)+cu),2));

% objective = objective + 1/dims.N*sum(norms(Bx*(Sx*x+cx*ones(1,dims.N)),2) + norms(Bu*(Su*u+cu*ones(1,dims.N)),2) + SCP.lambda*norms(E*v,2));

% Parameters In and Solutions Out:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    parameters = {A,Bm,Bp,F,E,r,H0,Hf,l0,lf,xbar,ubar,pbar,c,eta};
    solutions = {objective,Sx*x+cx*ones(1,dims.N),Su*u+cu*ones(1,dims.N),Sp*p+cp,E*v,vic,vtc};
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutiobs, {vs}];
    end
    
end
if SCP.FreeTime && SCP.AdaptiveMesh
    parameters = {A,Bm,Bp,H,E,r,H0,Hf,l0,lf,xbar,ubar,hbar,c,eta};
    solutions = {objective,Sx*x+cx*ones(1,dims.N),Su*u+cu*ones(1,dims.N),Sh*h+ch*ones(1,dims.N-1),E*v,vic,vtc};
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutiobs, {vs}];
    end
end
if ~SCP.FreeTime
    parameters = {A,Bm,Bp,E,r,H0,Hf,l0,lf,xbar,ubar,c,eta};
    solutions = {objective,Sx*x+cx*ones(1,dims.N),Su*u+cu*ones(1,dims.N),E*v,vic,vtc};
    if dims.ns ~= 0
        parameters = [parameters, {C,D,rp}];
        solutions = [solutiobs, {vs}];
    end
end

Controller = optimizer(constraints,objective,sdpsettings('solver','mosek'),parameters,solutions);

save([fileparts(mfilename('fullpath')) '\Controllers\' SCP.ControllerName],'Controller');

end

function [LinearizedDynamics,LinearizedConstraints] = Linearize(SCP,System,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition)

[x,u,c] = Variables(System);
[f,slack] = Dynamics(System,x,u,c);
[s] = NonConvexConstraints(System,x,u,c);
[g0] = InitialCondition(System,x,c);
[gf] = FinalCondition(System,x,c);

% Dynamics Jacobians
if ~SCP.FreeTime && ~SCP.AdaptiveMesh
    LinearizedDynamics.f = System.T*matlabFunction(f,"Vars",{x,u,c});
    LinearizedDynamics.A = System.T*matlabFunction(Jacobian(f,x),"Vars",{x,u,c});
    LinearizedDynamics.B = System.T*matlabFunction(Jacobian(f,u),"Vars",{x,u,c});
elseif SCP.FreeTime && ~SCP.AdaptiveMesh
    syms p;
    LinearizedDynamics.f = matlabFunction(p*f,"Vars",{x,u,p,c});
    LinearizedDynamics.A = matlabFunction(p*Jacobian(f,x),"Vars",{x,u,p,c});
    LinearizedDynamics.B = matlabFunction(p*Jacobian(f,u),"Vars",{x,u,p,c});
    LinearizedDynamics.F = matlabFunction(f,"Vars",{x,u,p,c});
elseif SCP.AdaptiveMesh
    syms h;
    LinearizedDynamics.f = matlabFunction(h*f,"Vars",{x,u,h,c});
    LinearizedDynamics.A = matlabFunction(h*Jacobian(f,x),"Vars",{x,u,h,c});
    LinearizedDynamics.B = matlabFunction(h*Jacobian(f,u),"Vars",{x,u,h,c});
    LinearizedDynamics.H = matlabFunction(f,"Vars",{x,u,h,c});
else
    error('Invalid Free Time and Adaptive Mesh Parameters Specified');
end
LinearizedDynamics.E = matlabFunction(Jacobian(x,slack));

% Non Convex Constraints Jacobians
if ~size(s) == 0
    LinearizedConstraints.s = matlabFunction(s,"Vars",{x,u,c});
    LinearizedConstraints.C = matlabFunction(Jacobian(s,x),"Vars",{x,u,c});
    LinearizedConstraints.D = matlabFunction(Jacobian(s,u),"Vars",{x,u,c});
end

% Initial Condition Jacobians
LinearizedConstraints.g0 = matlabFunction(g0,"Vars",{x,c});
LinearizedConstraints.H0 = matlabFunction(Jacobian(g0,x),"Vars",{x,c});

% Final Condition Jacobians
LinearizedConstraints.gf = matlabFunction(gf,"Vars",{x,c});
LinearizedConstraints.Hf = matlabFunction(Jacobian(gf,x),"Vars",{x,c});

end

function fx = Jacobian(f,x)

fx = [];

for i = 1:length(x)
    fx = [fx diff(f,x(i))];
end

end

function [Sx,Su,cx,cu] = ScalingMatrices(System,Scaling)

[xmin,xmax,umin,umax] = Scaling(System);

Sx = diag(xmax-xmin);
cx = xmin;

Su = diag(umax-umin);
cu = umin;

end

function [dims] = Dimensions(System,SCP,InitialGuess,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling)

[x,u,c] = Variables(System);
[f,slack] = Dynamics(System,x,u,c);
[lx,lu,ltf,Bx,Bu] = Cost();
[s] = NonConvexConstraints(System,x,u,c);
[g0] = InitialCondition(System,x,c);
[gf] = FinalCondition(System,x,c);
[xmin,xmax,umin,umax] = Scaling(System);

dims.nx = length(x);
if length(f) ~= dims.nx || length(lx) ~= dims.nx || size(Bx,2) ~= dims.nx || size(InitialGuess.x,1) ~= dims.nx || length(xmin) ~= dims.nx || length(xmax) ~= dims.nx
    error('Incommensurate Dimensions specified for state');
end

dims.nu = length(u);
if length(lu) ~= dims.nu || size(Bu,2) ~= dims.nu || size(InitialGuess.u,1) ~= dims.nu || length(umin) ~= dims.nu || length(umax) ~= dims.nu
    error('Incommensurate Dimensions specified for Input');
end

dims.N = size(InitialGuess.x,2);
if size(InitialGuess.u,2) ~= dims.N
    error('Incommensurate Dimensions specified for Nodes');
end

dims.nc = length(c);

dims.nic = length(g0);
dims.ntc = length(gf);
dims.nv = length(slack);
dims.ns = length(s);

if SCP.FreeTime && ~SCP.AdaptiveMesh
    if ~isfield(InitialGuess,'p')
        error('No Initial Time Guess Specified for Free Final Time Problem');
    end
end
if SCP.AdaptiveMesh
    if ~isfield(InitialGuess,'h')
        error('No Initial Mesh Guess Specified for Problem');
    end
    if length(InitialGuess.h) ~= dims.N-1
        error('Incommensurate Dimensions specified for initial mesh');
    end
end

end