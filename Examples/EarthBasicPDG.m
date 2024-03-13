function EarthBasicPDG

% Add SCP Path:
addpath('../Sequential Convex Programming/SCvx');
addpath('../Sequential Convex Programming/SCP Functions');
% Add Setup path:
addpath('../Setup');
% Add Plotting path:
addpath('../Plotting');

% Setup System:
System = EarthRocket;

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
SCP.alphah = 1;
% Sub Integration Values:
SCP.Nsub = 20; 
% Convergence Tolerance:
SCP.Tolerance = .01;
SCP.InfeasibilityTolerance = 10^-3;
SCP.MaximumIterationCount = 50;
% Time:
SCP.FreeTime = true;
SCP.AdaptiveMesh = false;
% Optimizer Saving:
SCP.RegenerateOptimizer = true;
SCP.OptimizerName = 'SixDoFPDG';
% Discretization:
SCP.Discretization = 'FOH';

% Boundary Conditions:
System.r0 = [0;200;300]; % Initial Position [m]
System.v0 = [-40;0;-40]; % Initial Velocity [m/s]
System.q0 = eul2quat([0,0,0])'; % Initial Orientation
System.w0 = [0;0;0]; % Initial Rotational Rate [rad/s]
System.rf = [0;0;0]; % Initial Position [m]
System.vf = [0;0;-5]; % Initial Velocity [m/s]
System.qf = eul2quat([0,0,0])'; % Initial Orientation
System.wf = [0;0;0]; % Initial Rotational Rate [rad/s]

% Scale:
rs = norm(System.r0);
ms = System.WetMass;
System.r0 = System.r0/rs;
System.v0 = System.v0/rs;
System.rf = System.rf/rs;
System.vf = System.vf/rs;
System.g = System.g/rs;
System.WetMass = System.WetMass/ms;
System.DryMass = System.DryMass/ms;
System.InertiaTensor = System.InertiaTensor/(ms*rs^2);
System.MaximumThrust = System.MaximumThrust/(ms*rs);
System.MinimumThrust = System.MinimumThrust/(ms*rs);
System.alpha = System.alpha*rs;
System.GimbaltoCOMDistance = System.GimbaltoCOMDistance/rs;
System.MaximumVelocity = System.MaximumVelocity/rs;

% Additional System Variables:
% Max Time:
t = 15; % Initial Time
System.tmax = inf;
System.tmin = 0;
% Boundary Conditions:
System.x0 = [System.r0;System.v0;System.q0;System.w0;System.WetMass];
System.xf = [System.rf;System.vf;System.qf;System.wf];

% Initial Guess:
N = 30; % Nodes
InitialGuess.x = System.x0*linspace(1,0,N) + [System.xf;System.DryMass]*linspace(0,1,N);
InitialGuess.u = (System.MaximumThrust-System.MinimumThrust)/2*[0;0;1;1]*ones(1,N);
InitialGuess.p = t;
InitialGuess.h = zeros(0,N);

% Call Sequential Convex Program
[Solution] = SCvx(SCP,System,InitialGuess,@Variables,@Cost,@Dynamics,@Convexities,@NonConvexConstraints,@InitialCondition,@FinalCondition,@Scaling,@Continuation);

for i = 1:length(Solution)
    x(:,:,i) = Solution{i}.x;
    u(:,:,i) = Solution{i}.u;
    p(i) = Solution{i}.p;
end

% Plot Results
x(1:6,:,:) = x(1:6,:,:)*rs;
x(14,:,:) = x(14,:,:)*ms;
PlotRocket(System,x(:,:,end),u(:,:,end));

end

function [x,u,p,h,c] = Variables(System)

% Create Symbolic State Variable:
syms rx ry rz vx vy vz q0 q1 q2 q3 w1 w2 w3 m;
x = [rx;ry;rz;vx;vy;vz;q0;q1;q2;q3;w1;w2;w3;m];

% Create Symbolic Input Variable:
syms u1 u2 u3 u4;
u = [u1;u2;u3;u4];

% Create Symbolic Parameter Variable:
syms t;
p = [t];

% Create Symbolic Time Varying Parameter Variable:
h = [];

% Create Symbolic Continuation Variable:
c = [];

end

%% Cost
% This Function gives the cost function of the optimization problem
% Total Cost: J = lx'*x(1) + lu'*u(1) + lp'*p + trapz(||Bx*x|| + ||Bu*u||)
%  Input:
%    System - 
%  Output:
%    phi - Terminal Cost
%    L - Running Cost
function [phi,L] = Cost(System,x,u,p,h,c)

% Terminal Cost:
phi = -100*x(14); % Minimize Fuel

% Running Cost:
L = 0; 

end

function [f,g,slack] = Dynamics(System,x,u,p,h,c)

% Rotation Matrices:
TB2L = @(q) [1-2*(q(3)^2+q(4)^2) 2*(q(2)*q(3)+q(4)*q(1)) 2*(q(2)*q(4)-q(3)*q(1));
             2*(q(2)*q(3)-q(4)*q(1)) 1-2*(q(2)^2+q(4)^2) 2*(q(3)*q(4)+q(2)*q(1));
             2*(q(2)*q(4)+q(3)*q(1)) 2*(q(3)*q(4)-q(2)*q(1)) 1-2*(q(2)^2+q(3)^2)]';

Omega = @(w) [0 -w(1) -w(2) -w(3);
              w(1) 0 w(3) -w(2);
              w(2) -w(3) 0 w(1);
              w(3) w(2) -w(1) 0];

% State Vectors:
r = x(1:3);
v = x(4:6);
q = x(7:10);
w = x(11:13);
m = x(14);

% Force, Moment:
F = u(1:3);
M = cross([0;0;-System.GimbaltoCOMDistance],F);

% Dynamics:
f = p(1)*[v;
         TB2L(q)*F/m + [0;0;-System.g];
         1/2*Omega(w)*q;
         System.InertiaTensor^-1*M-cross(w,w);%System.InertiaTensor^-1*(cross(w,System.InertiaTensor*w));
         -System.alpha*u(4)];

% Disturbance
g = zeros(14,0);

% Define Variables that can be slackened
slack = [x(1:13)];

end

function [g0] = InitialCondition(System,x0,u0,p,h0,c)

g0 = x0(1:14) - System.x0;

end

function [gf] = FinalCondition(System,xf,uf,p,hf,c)

gf = xf(1:13) - System.xf;

end

function [Constraints,Objective,Vars] = Convexities(System,x,u,p,h,c,Vars)

Constraints = [];

for k = 1:size(x,2)
    % State Vectors:
    r = x(1:3,k);
    v = x(4:6,k);
    q = x(7:10,k);
    w = x(11:13,k);
    m = x(14,k);
    
    % Nonconvex Inequality Constraint (s <= 0):
    Constraints = [Constraints;
                   m >= System.DryMass; % Mass Constraint
                   norm(r,2) <= tand(System.GlidescopeAngle)*r(3); % Glidescope Constraint
                   norm(w,2) <= System.MaximumRotationalRate; % Rotational Rate Constraint
                   norm(v,2) <= System.MaximumVelocity; % Velocity Constraint
                   norm(u(1:2,k),2) <= tand(System.MaximumGimbalAngle)*u(3,k); % Gimbal Constraint
                   norm(q(2:3),2) <= sqrt((1-cosd(System.MaximumPointingAngle))/2); % Pointing Angle Constraint
                   norm(u(1:3,k),2) <= u(4,k); % Slack Thrust Constraint
                   u(4,k) <= System.MaximumThrust; % Maximum Thrust Constraint
                   System.MinimumThrust <= u(4,k); % Minimum Thrust Constraint
                   ];
end

Objective = 0;

end

function [s] = NonConvexConstraints(System,x,u,p,h,c)

% Nonconvex Inequality Constraint (s <= 0):
s = [System.MinimumThrust^2 - (u(1)^2+u(2)^2+u(3)^2)];

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
function [xmin,xmax,umin,umax,pmin,pmax,hmin,hmax] = Scaling(System)

% State Scaling Limits
xmin = [min(System.r0,System.rf);
        -System.MaximumVelocity*ones(3,1);
        -ones(4,1);
        -System.MaximumRotationalRate*ones(3,1);
        System.DryMass];

xmax = [max(System.r0,System.rf);
        System.MaximumVelocity*ones(3,1);
        ones(4,1);
        System.MaximumRotationalRate*ones(3,1);
        System.WetMass];

% Input Scaling Limits
umin = [-System.MaximumThrust;
        -System.MaximumThrust;
        -System.MaximumThrust;
        -System.MaximumThrust];

umax = [System.MaximumThrust;
        System.MaximumThrust;
        System.MaximumThrust;
        System.MaximumThrust];

% Parameter Scaling Limits
pmin = [0];
pmax = [1];

% Time-varying Parameter Scaling Limits
hmin = zeros(0,1);
hmax = zeros(0,1);

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

