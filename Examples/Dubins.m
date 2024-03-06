function Dubins

% Add SCP Path:
addpath('../Numerical Methods/Sequential Convex Programming/SCvx');
addpath('../Numerical Methods/Sequential Convex Programming/SCP Functions');
% Plotting:
addpath('../Plotting');

% Setup System:
System.MaximumVelocity = 1;
System.MaximumRotationalRate = pi/6;
System.MaximumDistance = 10;
System.MinimumDistance = -10;
System.CarRadius = .5;

% Obstacles:
System.Obstacles(1).p = [5;4];
System.Obstacles(1).r = 3;
System.Obstacles(2).p = [-5;-4];
System.Obstacles(2).r = 3;
System.Obstacles(3).p = [0;0];
System.Obstacles(3).r = 2;

% Sequential Convex Program Variables:
% Trust Region Update Values:
SCP.betash = 2; 
SCP.betagr = 3.2;
% Accuracy Metric Regions:
SCP.rho0 = 0.0; % Case 1: rho < rho0
SCP.rho1 = .25; % Case 2: rho0 < rho < rho1
SCP.rho2 = .9; % Case 3: rho1 < rho < rh02
% Slack Penalization:
SCP.lambda = 10^5; % Dynamics Slack Penalty
SCP.lambdas = 10^5; % Non-Convexity Slack Penalty
SCP.lambdaic = 10^5; % Initial Condition Slack Penalty
SCP.lambdatc = 10^5; % Final Condition Slack Penalty
% Trust Region:
SCP.eta = 5; % Initial Trust Region
SCP.eta0 = 0; % Minimum Trust Region
SCP.eta1 = 100; % Maximum Trust Region
SCP.TrustRegionNorm = 1; % Norm to evaluate trust region
SCP.alphax = 1;
SCP.alphau = 1;
SCP.alphap = 1;
% Sub Integration Values:
SCP.Nsub = 20; 
% Convergence Tolerance:
SCP.Tolerance = 1;
SCP.InfeasibilityTolerance = 10^-3;
SCP.MaximumIterationCount = 50;
% Time:
SCP.FreeTime = true;
SCP.AdaptiveMesh = false;
% Optimizer Saving:
SCP.RegenerateOptimizer = true;
SCP.OptimizerName = 'DubinsCar';

% Boundary Conditions:
System.x0 = [-8;-8;0];
System.xf = [8;8;0];

% Non-dimensionalize:
rs = System.MaximumDistance-System.MinimumDistance;
System.MaximumVelocity = System.MaximumVelocity/rs;
System.MaximumDistance = System.MaximumDistance/rs;
System.MinimumDistance = System.MinimumDistance/rs;
System.CarRadius = System.CarRadius/rs;
System.x0(1:2) = System.x0(1:2)/rs;
System.xf(1:2) = System.xf(1:2)/rs;
System.Obstacles(1).p = System.Obstacles(1).p/rs;
System.Obstacles(1).r = System.Obstacles(1).r/rs;
System.Obstacles(2).p = System.Obstacles(2).p/rs;
System.Obstacles(2).r = System.Obstacles(2).r/rs;
System.Obstacles(3).p = System.Obstacles(3).p/rs;
System.Obstacles(3).r = System.Obstacles(3).r/rs;

% Additional System Variables:
% Max Time:
t = 30; % Initial Time
System.tmax = t*2;
System.tmin = t*.5;

% Initial Guess:
N = 30; % Nodes
InitialGuess.x = System.x0*linspace(1,0,N) + System.xf*linspace(0,1,N);
InitialGuess.u = zeros(2,N);
InitialGuess.p = t;

% Call Sequential Convex Program
[Solution] = SCvx(SCP,System,InitialGuess,@Variables,@Cost,@Dynamics,@ConvexConstraints,@NonConvexConstraints,@InitialCondition,@FinalCondition,@Scaling,@Continuation);

for i = 1:length(Solution)
    x(:,:,i) = Solution{i}.x;
    u(:,:,i) = Solution{i}.u;
    p(i) = Solution{i}.p;
end

x(1:2,:,:) = x(1:2,:,:)*rs;
u(1,:,i) = u(1,:,i)*rs;

% Plot Results
System.CarRadius = System.CarRadius*rs;
System.MaximumDistance = System.MaximumDistance*rs;
System.MinimumDistance = System.MinimumDistance*rs;
Boundary = [System.MinimumDistance System.MinimumDistance System.MaximumDistance System.MaximumDistance;
            System.MinimumDistance System.MaximumDistance System.MaximumDistance System.MinimumDistance];
PlotDubins(System,Boundary,x(:,:,end))

end

function [x,u,c] = Variables(System)

% Create Symbolic State Variable:
syms x y theta;
x = [x;y;theta];

% Create Symbolic Input Variable:
syms v w;
u = [v;w];

% Create Symbolic Continuation Variable:
c = sym([]);

end

%% Cost
% This Function gives the cost function of the optimization problem
% Total Cost: J = lx'*x(1) + lu'*u(1) + lp'*p + trapz(||Bx*x|| + ||Bu*u||)
%  Input:
%    System - 
%  Output:
%    phi - Terminal Cost
%    L - Running Cost
function [phi,L] = Cost(System,t,x,u,c)

% Terminal Cost:
phi = 10*t;

% Running Cost:
L = sym(0); 

end

function [f,slack] = Dynamics(System,x,u,c)


% Dynamics:
f = [u(1) * cos(x(3));
     u(1) * sin(x(3));
     u(2)];

% Define Variables that can be slackened
slack = [x(1:3)];

end

function [g0] = InitialCondition(System,x0,c)

g0 = x0(1:3) - System.x0;

end

function [gf] = FinalCondition(System,xf,c)

gf = xf(1:3) - System.xf;

end

function [Constraints] = ConvexConstraints(System,x,u,c)

% Nonconvex Inequality Constraint (s <= 0):
Constraints = [0 <= u(1);
               u(1) <= System.MaximumVelocity;
               u(2) <= System.MaximumRotationalRate;
               -System.MaximumRotationalRate <= u(2);
               x(1) <= System.MaximumDistance - System.CarRadius;
               System.MinimumDistance + System.CarRadius <= x(1);
               x(2) <= System.MaximumDistance - System.CarRadius;
               System.MinimumDistance + System.CarRadius <= x(2);
               ];

end

function [s] = NonConvexConstraints(System,x,u,c)

% Nonconvex Inequality Constraint (s <= 0):
s = [];
for i = 1:length(System.Obstacles)
    s = [s;
         System.CarRadius + System.Obstacles(i).r - norm(x(1:2)-System.Obstacles(i).p)];
end

end

%% Scaling
% This Function gives the minimum and maximum values of each solution
% variable so that the corresponding variables of the convex solver are
% approximately between 0 and 1.
%  Input:
%  Output:
%    xmin - Minimum State Vector
%    xmax - Maximum State Vector
%    umin - Minimum Input Vector
%    umax - Maximum Input Vector
%    pmin - Minimum Parameter Vector
%    pmax - Maximum Parameter Vector
function [xmin,xmax,umin,umax] = Scaling(System)

% State Scaling Limits
xmin = [System.MinimumDistance;
        System.MinimumDistance;
        -2*pi];

xmax = [System.MaximumDistance;
        System.MaximumDistance;
        2*pi];

% Input Scaling Limits
umin = [0;
        -System.MaximumRotationalRate];

umax = [System.MaximumVelocity;
        System.MaximumRotationalRate];

end

%% Continuation
% This Function gives applies continuation to the system parameters that
% govern the trajectory optimization problem. A common continuation method
% concerns the maximum thrust of a spacecraft
%  Input:
%    System - System variable structure
%    gamma - Continuation parameter [0,1]
%  Output:
%    c - Continuation Parameter Value
function [c] = Continuation(System,gamma)

% Apply Continuation to Thrust:
c = [];

end

