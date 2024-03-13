function ApolloPDG

% Add SCP Path:
addpath('../Sequential Convex Programming/SCvx');
addpath('../Sequential Convex Programming/SCP Functions');
% Add Setup path:
addpath('../Setup');
% Add Plotting path:
addpath('../Plotting');

% Setup System:
System = ApolloDescentLander;

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
SCP.Nsub = 50; 
% Convergence Tolerance:
SCP.Tolerance = .01;
SCP.InfeasibilityTolerance = 10^-3;
SCP.MaximumIterationCount = 50;
% Time:
SCP.FreeTime = true;
SCP.AdaptiveMesh = false;
% Optimizer Saving:
SCP.RegenerateOptimizer = true;
SCP.OptimizerName = 'ApolloPDG';
% Discretization:
SCP.Discretization = 'FOH';

% Boundary Conditions:
System.r0 = [-200;0;300]; % Initial Position [m]
System.v0 = [20;0;0]; % Initial Velocity [m/s]
System.q0 = eul2quat([0,-pi/2,0])'; % Initial Orientation
System.w0 = [0;0;0]; % Initial Rotational Rate [rad/s]
System.rf = [0;0;0]; % Initial Position [m]
System.vf = [0;0;0]; % Initial Velocity [m/s]
System.qf = eul2quat([0,0,0])'; % Initial Orientation
System.wf = [0;0;0]; % Initial Rotational Rate [rad/s]

% Additional System Variables:
t = 60; % Final Time
System.tmax = inf;
System.tmin = 0;

% Scale:
rs = norm(System.r0);
ms = System.WetMass;
ts = t;
t = t/ts;
System.tmax = System.tmax/ts;
System.tmin = System.tmin/ts;
System.r0 = System.r0/rs;
System.v0 = System.v0*ts/rs;
System.rf = System.rf/rs;
System.vf = System.vf*ts/rs;
System.g = System.g*ts^2/rs;
System.WetMass = System.WetMass/ms;
System.DryMass = System.DryMass/ms;
System.InertiaTensor = System.InertiaTensor/(ms*rs^2);
System.DPSMaximumThrust = System.DPSMaximumThrust*ts^2/(ms*rs);
System.DPSMinimumThrust = System.DPSMinimumThrust*ts^2/(ms*rs);
System.RCSMaximumThrust = System.RCSMaximumThrust*ts^2/(ms*rs);
System.alphaDPS = System.alphaDPS*rs/ts;
System.alphaRCS = System.alphaRCS*rs/ts;
System.MaximumVelocity = System.MaximumVelocity*ts/rs;
System.DPSPropellantMass = System.DPSPropellantMass/ms;
System.RCSPropellantMass = System.RCSPropellantMass/ms;
for i = 1:length(System.RCSThruster)
    System.RCSThruster(i).Location = System.RCSThruster(i).Location/rs;
end

% Boundary Conditions:
System.x0 = [System.r0;System.v0;System.q0;System.w0;System.DPSPropellantMass;System.RCSPropellantMass];
%PlotApollo(System,System.x0,ones(13,1));
System.xf = [System.rf;System.vf;System.qf;System.wf];

% Initial Guess:
N = 30; % Nodes
InitialGuess.x = System.x0*linspace(1,0,N) + [System.xf;System.DPSPropellantMass;System.RCSPropellantMass]*linspace(0,1,N);
InitialGuess.u = (System.DPSMaximumThrust+System.DPSMinimumThrust)/2*[1;zeros(12,1)]*ones(1,N);
InitialGuess.p = t;
InitialGuess.h = zeros(0,N);

% Call Sequential Convex Program
[Solution] = SCvx(SCP,System,InitialGuess,@Variables,@Cost,@Dynamics,@Convexities,@NonConvexConstraints,@InitialCondition,@FinalCondition,@Scaling,@Continuation);

for i = 1:length(Solution)
    x(:,:,i) = Solution{i}.x;
    u(:,:,i) = Solution{i}.u;
    p(i) = Solution{i}.p;
end
% for i = 1:length(Solution)-1
%     J(i) = Solution{i}.J;
%     L(i) = Solution{i+1}.J;
% end

% Plot Results
x(1:6,:,:) = x(1:6,:,:)*rs;
x(14:15,:,:) = x(14:15,:,:)*ms;
u(:,:,:) = u(:,:,:)*ms*rs/ts^2;
System.MaximumThrust = System.DPSMaximumThrust*ms*rs/ts^2;
PlotApollo(System,x(:,:,end),u(:,:,end));

end

function [x,u,p,h,c] = Variables(System)

% Create Symbolic State Variable:
syms rx ry rz vx vy vz q0 q1 q2 q3 w1 w2 w3 mDPS mRCS;
x = [rx;ry;rz;vx;vy;vz;q0;q1;q2;q3;w1;w2;w3;mDPS;mRCS];

% Create Symbolic Input Variable:
syms uDPS uRCS1 uRCS2 uRCS3 uRCS4 uRCS5 uRCS6 uRCS7 uRCS8 uRCS9 uRCS10 uRCS11 uRCS12;
u = [uDPS;uRCS1;uRCS2;uRCS3;uRCS4;uRCS5;uRCS6;uRCS7;uRCS8;uRCS9;uRCS10;uRCS11;uRCS12];

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

% Force, Moment:
F = [0;0;1]*u(1);
M = 0;
for i = 1:length(System.RCSThruster)
    F = F - System.RCSThruster(i).Direction*u(i+1);
    M = M - cross(System.RCSThruster(i).Location,System.RCSThruster(i).Direction)*u(i+1);
end

% Dynamics:
f = p(1)*[v;
         TB2L(q)*F/(System.DryMass+x(14)+x(15)) + [0;0;-System.g];
         1/2*Omega(w)*q;
         System.InertiaTensor^-1*M;%-System.InertiaTensor^-1*(cross(w,System.InertiaTensor*w));
         -System.alphaDPS*u(1);
         -System.alphaRCS*sum(u(2:13))];

% Disturbance
g = zeros(15,0);

% Define Variables that can be slackened
slack = [x(1:15)];

end

function [g0] = InitialCondition(System,x0,u0,p,h0,c)

g0 = x0(1:15) - System.x0;

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
    mDPS = x(14,k);
    mRCS = x(15,k);
    
    % Nonconvex Inequality Constraint (s <= 0):
    Constraints = [Constraints;
                   mRCS >= 0; % Mass Constraint
                   mDPS >= 0; 
                   %norm(r,2) <= tand(System.GlidescopeAngle)*r(3); % Glidescope Constraint
                   %norm(w,2) <= System.MaximumRotationalRate; % Rotational Rate Constraint
                   %norm(v,2) <= System.MaximumVelocity; % Velocity Constraint
                   norm(q(2:3),2) <= sqrt((1-cosd(System.MaximumPointingAngle))/2); % Pointing Angle Constraint
                   u(1,k) >= System.DPSMinimumThrust;
                   u(1,k) <= System.DPSMaximumThrust
                   r(3) >= 0];

    for j = 2:13
        Constraints = [Constraints;
                       u(j,k) >= 0;
                       u(j,k) <= System.RCSMaximumThrust];
    end

end

Objective = 0;

end

function [s] = NonConvexConstraints(System,x,u,p,h,c)

% Nonconvex Inequality Constraint (s <= 0):
s = [];

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
        0;
        0];

xmax = [max(System.r0,System.rf);
        System.MaximumVelocity*ones(3,1);
        ones(4,1);
        System.MaximumRotationalRate*ones(3,1);
        System.DPSPropellantMass;
        System.RCSPropellantMass];

% Input Scaling Limits
umin = [-System.DPSMaximumThrust;
        -System.RCSMaximumThrust*ones(12,1)];

umax = [System.DPSMaximumThrust;
        System.RCSMaximumThrust*ones(12,1)];

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

