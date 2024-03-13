
function [L,Solution] = SCvx(SCP,System,InitialGuess,Variables,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition,Scaling,Continuation)

% Algorithm Start
tic;

% Determine Problem Dimensions
[dims] = Dimensions(System,SCP,InitialGuess,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling);

% Initialize Continuation Parameter:
if dims.nc == 0
    gamma = 1;
    InitialGuess.c = [];
else
    gamma = 0;
    InitialGuess.c = Continuation(System,gamma);
end

% Scaling Matrices:
[ScalingMatrices] = Scale(SCP,System,Scaling);

% Start Iteration:
disp(['----------- Iteration ' num2str(1) ' -----------']);

% Find Jacobians
fprintf('Linearizing Problem...');
[LinearizedDynamics,LinearizedConstraints] = Linearize(SCP,System,Variables,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition);
fprintf('Done\n');

% Initialize MOSEK:
[~, MosekCodes] = mosekopt('symbcon echo(0)');

% Discretize
fprintf('Discretizing Problem...');
[DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,InitialGuess,LinearizedDynamics,LinearizedConstraints,dims);
fprintf('Done\n');

% Parse the Problem
fprintf('Parsing Problem...');
[ConvexProblem] = Parse(SCP,System,InitialGuess,Cost,DiscreteDynamics,DiscreteConstraints,ScalingMatrices,dims,MosekCodes);
fprintf('Done\n');

% Nonlinear Cost:
Jbar = NonLinearCost(SCP,InitialGuess,LinearizedConstraints,Cost,delta,dims);

Solution{1} = InitialGuess;
Solution{1}.Tolerance = NaN;
Solution{1}.delta = delta;
iter = 2;

while true
    
    % Solve Problem as a SOCP
    fprintf('Solving Convex Subproblem..');
    [Lstar,Solution{iter}] = ConvexSolve(SCP,ConvexProblem,ScalingMatrices,Cost,dims);
    Solution{iter}.c = Solution{iter-1}.c;
    fprintf('Done\n');

    % Check Convergence Criterion
    if iter ~=1
        Solution{iter}.Tolerance = abs(Jbar-Lstar);
        Solution{iter}.InfeasibilityTolerance = max(vecnorm(delta));
        if Solution{iter}.Tolerance <= SCP.Tolerance && gamma == 1 && Solution{iter}.InfeasibilityTolerance <= SCP.InfeasibilityTolerance
            disp(['SCP Algorithm Converged in ' num2str(toc) ' seconds']);
            break;
        else
            disp(['Current Iteration: ' num2str(iter) ' Current Tolerance: ' num2str(abs(Jbar-Lstar)) ' Current Trust Region: ' num2str(SCP.eta)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
        end
    
        if iter >= SCP.MaximumIterationCount
            disp(['Maximum Iteration Count Reached. Current Tolerance: ' num2str(abs(Jbar-Lstar))]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
            break
        end
    end
    
    % Redescretize with solution
    fprintf('Discretizing Problem...');
    [DiscreteDynamicsStar,DiscreteConstraintsStar,delta] = Discretize(SCP,NewReference,LinearizedDynamics,LinearizedConstraints,dims);
    fprintf('Done\n');

    % New Nonlinear Cost
    Jstar = NonLinearCost(SCP,NewReference,LinearizedConstraints,Cost,delta,dims);
    
    % Calculate Convexification Accuracy
    rho = (Jbar - Jstar)/(Jbar-Lstar);
    disp(['Current rho: ' num2str(rho)]);

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

    % Parse the Problem
    fprintf('Parsing Problem...');
    [ConvexProblem] = Parse(SCP,System,InitialGuess,Cost,DiscreteDynamics,DiscreteConstraints,ScalingMatrices,dims,MosekCodes);
    fprintf('Done\n');

    iter = iter + 1;
    
end


end

function J = NonLinearCost(SCP,Solution,LinearizedConstraints,Cost,delta,dims)

% Cost:
[lx,lu,ltf,Bx,Bu] = Cost();

% Terminal Cost:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    phi = lx'*Solution.x(:,end) + lu'*Solution.u(:,end) + ltf*Solution.p + SCP.lambdaic*norm(LinearizedConstraints.g0(Solution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(Solution.x(:,end)),2);
end
if SCP.FreeTime && SCP.AdaptiveMesh
    phi = lx'*Solution.x(:,end) + lu'*Solution.u(:,end) + ltf*sum(Solution.h) + SCP.lambdaic*norm(LinearizedConstraints.g0(Solution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(Solution.x(:,end)),2);
end
if ~SCP.FreeTime
    phi = lx'*Solution.x(:,end) + lu'*Solution.u(:,end) + SCP.lambdaic*norm(LinearizedConstraints.g0(Solution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(Solution.x(:,end)),2);
end

% Running Cost:
if dims.ns ~= 0
    Gamma = 1/dims.N*(sum(vecnorm(Bx*Solution.x)) + sum(vecnorm(Bu*Solution.u)) + SCP.lambda*sum(vecnorm(delta)) + SCP.lambdas*vecnorm(max(LinearizedConstraints.s(ReferenceSolution.x,ReferenceSolution.u),0)));
else
    Gamma = 1/dims.N*(sum(vecnorm(Bx*Solution.x)) + sum(vecnorm(Bu*Solution.u)) + SCP.lambda*sum(vecnorm(delta)));
end

% Total Cost:
J = phi + Gamma;

end

function [Solution] = DecomposeResult(Result)

% Extract Primal Solution
z = Result.sol.itr.xx;

% Decompose Solution
n = 1;
Solution.x = ScalingMatrices.Sx*reshape(z(n:n+dims.nx*dims.N-1),dims.nx,dims.N)+ScalingMatrices.cx; n = n+dims.nx*dims.N; % State
Solution.u = ScalingMatrices.Su*reshape(z(n:n+dims.nu*dims.N-1),dims.nu,dims.N)+ScalingMatrices.cu; n = n+dims.nu*dims.N; % Input
if SCP.AdaptiveMesh
    Solution.h = ScalingMatrices.Sh*z(n:n+dims.N-1-1)+ScalingMatrices.ch; n = n+dims.N-1; % Mesh Parameter
end
if SCP.FreeTime && ~SCP.AdaptiveMesh
    Solution.p = ScalingMatrices.Sp*z(n)+ScalingMatrices.cp; n = n+1; % Time Parameter
end
Solution.v = reshape(z(n:n+dims.nv*(dims.N-1)-1),[dims.nv,dims.N-1]); n = n+dims.nv*(dims.N-1); % Dynamic Infeasibility
Solution.vs = reshape(z(n:n+dims.ns*dims.N-1),[dims.ns,dims.N]); n = n+dims.ns*dims.N; % Nonconvex constraint Infeasibility
Solution.vic = z(n:n+dims.nic-1); n = n+dims.nic; % Initial Condition Dynamic Infeasibility
Solution.vtc = z(n:n+dims.ntc-1); n = n+dims.ntc; % Final Condition Dynamic Infeasibility
Solution.mux = z(n:n+dims.N-1); n = n+dims.N; % State Trust Region Constraint
Solution.muu = z(n:n+dims.N-1); n = n+dims.N; % Input Trust Region Constraint
if SCP.AdaptiveMesh
    Solution.muh = z(n:n+dims.N-1-1); n = n+dims.N-1; % Mesh Trust Region Constraint
end
if SCP.FreeTime && ~SCP.AdaptiveMesh
    Solution.mup = z(n); n = n+1; % Time Trust Region Constraint
end
[lx,lu,~,Bx,Bu] = Cost();
if ~all(Bx(:)==0)
    Solution.sigmax = z(n:n+dims.N-1); n = n+dims.N; % State Cost Constraint
end
if ~all(Bu(:)==0)
    Solution.sigmau = z(n:n+dims.N-1); n = n+dims.N; % State Cost Constraint
end
Solution.sigmav = z(n:n+dims.N-1-1); n = n+dims.N-1; % Dynamic Infeasibility Cost Constraint
if dims.ns ~= 0
    Solution.sigmavs = z(n:n+dims.N-1); n = n+dims.N; % Non Convex Constraint Infeasibility Cost Constraint
end
Solution.sigmaic = z(n); n = n+1; % Initial Condition Infeasibility Cost Constraint
Solution.sigmatc = z(n);  % Initial Condition Infeasibility Cost Constraint

L = Result.sol.itr.pobjval+ConvexProblem.cf; % Linear Cost

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