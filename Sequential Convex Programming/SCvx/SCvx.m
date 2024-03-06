
function [Solution] = SCvx(SCP,System,InitialGuess,Variables,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition,Scaling,Continuation)

% Algorithm Start
tic;

% Determine Problem Dimensions
[dims] = Dimensions(System,SCP,InitialGuess,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling);

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

% Parse the Problem
fprintf('Parsing Problem...');
[ConvexProblemOptimizer] = YalmipParseSCvx(SCP,System,ConvexConstraints,ScalingMatrices,dims);
fprintf('Done\n');

% Start Iteration:
disp(['----------- Iteration ' num2str(1) ' -----------']);

% Find Jacobians
fprintf('Linearizing Problem...');
[LinearizedCost,LinearizedDynamics,LinearizedConstraints] = Linearize(SCP,System,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition);
fprintf('Done\n');

% Discretize
fprintf('Discretizing Problem...');
[DiscreteCost,DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,System,InitialGuess,LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims);
fprintf('Done\n');

% Initialize First Solution
Solution{1} = InitialGuess;
Solution{1}.Tolerance = NaN;
Solution{1}.delta = delta;

% Nonlinear Cost:
Solution{1}.J = NonLinearCost(SCP,InitialGuess,LinearizedCost,LinearizedConstraints,delta,dims);

SolutionIterate = 2;
ActualIterate = 2;

while true

    % Start Iteration:
    disp(['----------- Iteration ' num2str(SolutionIterate) ' -----------']);
    
    % Solve Problem
    fprintf('Solving Convex Subproblem..');
    [Solution{SolutionIterate}] = ConvexSolve(SCP,Solution{SolutionIterate-1},ConvexProblemOptimizer,DiscreteCost,DiscreteDynamics,DiscreteConstraints,dims);
    fprintf('Done\n');

    % Check Convergence Criterion
    Solution{SolutionIterate}.Tolerance = Solution{SolutionIterate-1}.J-Solution{SolutionIterate}.L;
    Solution{SolutionIterate}.InfeasibilityTolerance = max(vecnorm(delta));
    if abs(Solution{SolutionIterate}.Tolerance) <= SCP.Tolerance && gamma == 1 && Solution{SolutionIterate}.InfeasibilityTolerance <= SCP.InfeasibilityTolerance
        disp(['SCP Algorithm Converged in ' num2str(toc) ' seconds']);
        break;
    else
        disp(['Current Tolerance: ' num2str(Solution{SolutionIterate}.Tolerance) ' Current Trust Region: ' num2str(SCP.eta)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
    end

    if ActualIterate >= SCP.MaximumIterationCount
        disp(['Maximum Iteration Count Reached. Current Tolerance: ' num2str(Solution{SolutionIterate}.Tolerance)]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
        break
    end
    
    % Redescretize with solution
    fprintf('Discretizing Problem...');
    [NewDiscreteCost,NewDiscreteDynamics,NewDiscreteConstraints,delta] = Discretize(SCP,System,Solution{SolutionIterate},LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims);
    fprintf('Done\n');

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

% Final Time:
if SCP.FreeTime && ~SCP.AdaptiveMesh
    T = Solution.p;
elseif SCP.AdaptiveMesh
    T = sum(Solution.h);
else
    T = System.T;
end

% Terminal Cost:
% Basic Terminal Cost:
phi = LinearizedCost.phi(T,Solution.x(:,dims.N),Solution.c); 
% Augmented Terminal Cost:
phi = phi + SCP.lambdaic*norm(LinearizedConstraints.g0(Solution.x(:,1)),2) + SCP.lambdatc*norm(LinearizedConstraints.gf(Solution.x(:,dims.N)),2);

% Running Cost:
L = 0;
for k = 1:dims.N-1
    % Basic Running Cost:
    L = L + 1/dims.N*LinearizedCost.L(T,Solution.x(:,k),Solution.u(:,k),Solution.c);
    % Augmented Running Cost:
    L = L + 1/dims.N*(SCP.lambda*norm(delta(:,k),1));
    if dims.ns ~= 0
        L = L + 1/dims.N*(SCP.lambdas*sum(max(0,LinearizedConstraints.s(Solution.x(:,k),Solution.u(:,k))),1));
    end
end

% Total Cost:
J = phi + L;

end

function [dims] = Dimensions(System,SCP,InitialGuess,Variables,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition,Scaling)

[x,u,c] = Variables(System);
[f,slack] = Dynamics(System,x,u,c);
[s] = NonConvexConstraints(System,x,u,c);
[g0] = InitialCondition(System,x,c);
[gf] = FinalCondition(System,x,c);
[xmin,xmax,umin,umax] = Scaling(System);

dims.nx = length(x);
if length(f) ~= dims.nx || size(InitialGuess.x,1) ~= dims.nx || length(xmin) ~= dims.nx || length(xmax) ~= dims.nx
    error('Incommensurate Dimensions specified for state');
end

dims.nu = length(u);
if size(InitialGuess.u,1) ~= dims.nu || length(umin) ~= dims.nu || length(umax) ~= dims.nu
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