
function [Solution] = SCvx(SCP,System,InitialGuess,Variables,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition,Scaling,Continuation)

% Algorithm Start
tic;

% Determine Problem Dimensions
[dims] = Dimensions(System,InitialGuess,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling);

% Initialize Continuation Parameter:
if dims.nc == 0
    gamma = 1;
    InitialGuess.c = [];
else
    gamma = 0;
    InitialGuess.c = Continuation(System,gamma);
end

% Scaling Matrices:
[ScalingMatrices] = Scale(System,Scaling);

% Parse the Problem
fprintf('Parsing Problem...');
[ConvexProblemOptimizer,Fields] = YalmipParseSCvx(SCP,System,ConvexConstraints,ScalingMatrices,dims);
fprintf('Done\n');

% Find Jacobians
fprintf('Linearizing Problem...');
[LinearizedCost,LinearizedDynamics,LinearizedConstraints] = Linearize(System,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition);
fprintf('Done\n');

% Start Iteration:
disp(['----------- Iteration ' num2str(1) ' -----------']);

% Discretize
fprintf('Discretizing Problem...');
[DiscreteCost,DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,InitialGuess,LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims);
fprintf('Done\n');

% Initialize First Solution
Solution{1} = InitialGuess;
Solution{1}.Tolerance = NaN;
Solution{1}.delta = delta;
Solution{1}.Trajectory = DiscreteDynamics.Trajectory;

% Nonlinear Cost:
Solution{1}.J = NonLinearCost(SCP,InitialGuess,LinearizedCost,LinearizedConstraints,delta,dims);

SolutionIterate = 2;
ActualIterate = 2;

while true

    % Start Iteration:
    disp(['----------- Iteration ' num2str(SolutionIterate) ' -----------']);
    
    % Solve Problem
    fprintf('Solving Convex Subproblem..');
    [Solution{SolutionIterate}] = ConvexSolve(SCP,Solution{SolutionIterate-1},ConvexProblemOptimizer,DiscreteCost,DiscreteDynamics,DiscreteConstraints,Fields,dims);
    fprintf('Done\n');

    % Redescretize with solution
    fprintf('Discretizing Problem...');
    [NewDiscreteCost,NewDiscreteDynamics,NewDiscreteConstraints,delta] = Discretize(SCP,Solution{SolutionIterate},LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims);
    fprintf('Done\n');

    % Add solution Variables:
    Solution{SolutionIterate}.delta = delta;
    Solution{SolutionIterate}.Trajectory = DiscreteDynamics.Trajectory;

    % Check Convergence Criterion
    Solution{SolutionIterate}.Tolerance = abs(Solution{SolutionIterate-1}.J-Solution{SolutionIterate}.L);
    Solution{SolutionIterate}.InfeasibilityTolerance = max(vecnorm(delta));
    if Solution{SolutionIterate}.Tolerance <= SCP.Tolerance && gamma == 1 && Solution{SolutionIterate}.InfeasibilityTolerance <= SCP.InfeasibilityTolerance
        disp(['SCP Algorithm Converged in ' num2str(toc) ' seconds']);
        break;
    else
        disp(['Current Tolerance: ' num2str(Solution{SolutionIterate}.Tolerance) ' Current Trust Region: ' num2str(SCP.eta)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
    end

    if ActualIterate >= SCP.MaximumIterationCount
        disp(['Maximum Iteration Count Reached. Current Tolerance: ' num2str(Solution{SolutionIterate}.Tolerance)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
        break
    end

    % New Nonlinear Cost
    Solution{SolutionIterate}.J = NonLinearCost(SCP,Solution{SolutionIterate},LinearizedCost,LinearizedConstraints,delta,dims);
    
    % Calculate Convexification Accuracy
    Solution{SolutionIterate}.rho = (Solution{SolutionIterate-1}.J - Solution{SolutionIterate}.J)/(Solution{SolutionIterate-1}.J-Solution{SolutionIterate}.L);
    disp(['Current rho: ' num2str(Solution{SolutionIterate}.rho)]);

    % Update trust region
    if Solution{SolutionIterate}.rho < SCP.rho0
        % Update eta, dont update trajectory
        SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
    else
        if Solution{SolutionIterate}.rho < SCP.rho1
            % Update eta
            SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
        end
        % If rho1 < rho <= rho2, do not update eta
        if Solution{SolutionIterate}.rho >= SCP.rho2
            % Update eta
            SCP.eta = min(SCP.eta1,SCP.betagr*SCP.eta);
        end
        % Update Discretiszation
        DiscreteCost = NewDiscreteCost;
        DiscreteDynamics = NewDiscreteDynamics;
        DiscreteConstraints = NewDiscreteConstraints;
        % Move to Next Iterate:
        SolutionIterate = SolutionIterate + 1;
    end
    
    % Update Actual Iterate:
    ActualIterate = ActualIterate + 1;
    
end

end

function J = NonLinearCost(SCP,Solution,LinearizedCost,LinearizedConstraints,delta,dims)

% Terminal Cost:
% Basic Terminal Cost:
phi = LinearizedCost.phi(Solution.x(:,dims.N),Solution.u(:,dims.N),Solution.p,Solution.h(:,end),Solution.c); 
% Augmented Terminal Cost:
g0 = LinearizedConstraints.g0(Solution.x(:,1),Solution.u(:,1),Solution.p,Solution.h(:,1),Solution.c);
gf = LinearizedConstraints.gf(Solution.x(:,end),Solution.u(:,end),Solution.p,Solution.h(:,end),Solution.c);
phi = phi + SCP.lambdaic*norm(g0,2) + SCP.lambdatc*norm(gf,2);

% Running Cost:
L = 0;
for k = 1:dims.N
    % Basic Running Cost:
    L = L + 1/dims.N*LinearizedCost.L(Solution.x(:,k),Solution.u(:,k),Solution.p,Solution.h(:,k),Solution.c);
    % Augmented Running Cost:
    if k < dims.N
        L = L + 1/dims.N*(SCP.lambda*norm(delta(:,k),1));
    end
    if dims.ns ~= 0
        L = L + 1/dims.N*(SCP.lambdas*sum(max(0,LinearizedConstraints.s(Solution.x(:,k),Solution.u(:,k),Solution.p,Solution.h(:,k),Solution.c)),1));
    end
end

% Total Cost:
J = phi + L;

end

function [dims] = Dimensions(System,InitialGuess,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling)

[x,u,p,h,c] = Variables(System);
[f,g,slack] = Dynamics(System,x,u,p,h,c);
[s] = NonConvexConstraints(System,x,u,p,h,c);
[g0] = InitialCondition(System,x,u,p,h,c);
[gf] = FinalCondition(System,x,u,p,h,c);
[xmin,xmax,umin,umax,pmin,pmax,hmin,hmax] = Scaling(System);

dims.nx = length(x);
if length(f) ~= dims.nx || size(InitialGuess.x,1) ~= dims.nx || length(xmin) ~= dims.nx || length(xmax) ~= dims.nx
    error('Incommensurate Dimensions specified for state');
end

dims.nu = length(u);
if size(InitialGuess.u,1) ~= dims.nu || length(umin) ~= dims.nu || length(umax) ~= dims.nu
    error('Incommensurate Dimensions specified for Input');
end

dims.np = length(p);
if length(InitialGuess.p) ~= dims.np || length(pmin) ~= dims.np || length(pmax) ~= dims.np
    error('Incommensurate Dimensions Specified for Parameter');
end

dims.nh = length(h);
if size(InitialGuess.h,1) ~= dims.nh || length(hmin) ~= dims.nh || length(hmax) ~= dims.nh
    error('Incommensurate Dimensions Specified for Time-Varying Parameter');
end

dims.ng = size(g,2);

dims.N = size(InitialGuess.x,2);
if size(InitialGuess.x,2) ~= dims.N || size(InitialGuess.u,2) ~= dims.N || size(InitialGuess.h,2) ~= dims.N
    error('Incommensurate Dimensions specified for Nodes');
end

dims.nc = length(c);
dims.nic = length(g0);
dims.ntc = length(gf);
dims.nv = length(slack);
dims.ns = length(s);

end