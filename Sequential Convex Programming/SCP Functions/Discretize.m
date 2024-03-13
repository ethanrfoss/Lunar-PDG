function [DiscreteCost,DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,Reference,LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims)

% Discretize Dynamics:
[DiscreteDynamics,delta] = DiscretizeDynamics(SCP,Reference,LinearizedDynamics,dims);

% Discretize Nonconvex Constraints:
DiscreteConstraints = DiscretizeNonConvexConstraints(Reference,LinearizedConstraints,dims);

% Discretize Boundary Conditions:
DiscreteConstraints.H0 = LinearizedConstraints.H0(Reference.x(:,1),Reference.u(:,1),Reference.p,Reference.h(:,1),Reference.c);
DiscreteConstraints.I0 = LinearizedConstraints.I0(Reference.x(:,1),Reference.u(:,1),Reference.p,Reference.h(:,1),Reference.c);
DiscreteConstraints.J0 = LinearizedConstraints.J0(Reference.x(:,1),Reference.u(:,1),Reference.p,Reference.h(:,1),Reference.c);
DiscreteConstraints.K0 = LinearizedConstraints.K0(Reference.x(:,1),Reference.u(:,1),Reference.p,Reference.h(:,1),Reference.c);
DiscreteConstraints.l0 = LinearizedConstraints.g0(Reference.x(:,1),Reference.u(:,1),Reference.p,Reference.h(:,1),Reference.c) - DiscreteConstraints.H0*Reference.x(:,1) - DiscreteConstraints.I0*Reference.u(:,1) - DiscreteConstraints.J0*Reference.p - DiscreteConstraints.K0*Reference.h(:,1);

DiscreteConstraints.Hf = LinearizedConstraints.Hf(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteConstraints.If = LinearizedConstraints.If(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteConstraints.Jf = LinearizedConstraints.Jf(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteConstraints.Kf = LinearizedConstraints.Kf(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteConstraints.lf = LinearizedConstraints.gf(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c) - DiscreteConstraints.Hf*Reference.x(:,end) - DiscreteConstraints.If*Reference.u(:,end) - DiscreteConstraints.Jf*Reference.p - DiscreteConstraints.Kf*Reference.h(:,end);

% Discretize Cost:
[DiscreteCost] = DiscretizeCost(Reference,LinearizedCost,dims);

end

function [DiscreteDynamics,delta] = DiscretizeDynamics(SCP,Reference,LinearizedDynamics,dims)

% Preallocate:
delta = zeros(dims.nx,dims.N-1);
DiscreteDynamics.A = zeros(dims.nx,dims.nx,dims.N-1);
DiscreteDynamics.Bm = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.Bp = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.C = zeros(dims.nx,dims.np,dims.N-1);
DiscreteDynamics.D = zeros(dims.nx,dims.nh,dims.N-1);
DiscreteDynamics.G = zeros(dims.nx,dims.nx,dims.N-1);
DiscreteDynamics.Trajectory = [];

% Flattening Shape:
shape = [dims.nx,1;
         dims.nx,dims.nx;
         dims.nx,dims.nu;
         dims.nx,dims.nu;
         dims.nx,dims.np;
         dims.nx,dims.nh;
         dims.nx,1;
         dims.nx,dims.nx];

% Time:
t = linspace(0,1,dims.N);

% For Loop to Compute Discrete Time Dynamics at Each Time Step
for k = 1:dims.N-1
    
    % Initial Conditions
    phi0 = eye(dims.nx); % Initial Condition for A
    PBm0 = zeros(dims.nx,dims.nu); % Initial Condition for Bm
    PBp0 = zeros(dims.nx,dims.nu); % Initial Condition for Bp
    PC0 = zeros(dims.nx,dims.np);  % Initial Condition for C
    PD0 = zeros(dims.nx,dims.nh); % Initial Condition for D
    Pr0 = zeros(dims.nx,1); % Initial Condition for r
    PG0 = zeros(dims.nx,dims.nx); % Initial Condition for G

    % Vectorize Initial Conditions
    P0 = Flatten({Reference.x(:,k),phi0,PBm0,PBp0,PC0,PD0,Pr0,PG0},shape);
    
    % Integrate from current time step to next:
    tspan = linspace(t(k),t(k+1),SCP.Nsub);
    [P,Trajectory] = RK4(@(tc,P) Derivatives(SCP,P,tc,Reference.u(:,k+1),Reference.u(:,k),Reference.p,Reference.h(:,k),Reference.c,t(k+1),t(k),shape,LinearizedDynamics),tspan,P0);
    
    % Unvectorize Integrated Result
    M = Unflatten(P,shape);
    xf = M{1}; 
    Pphi = M{2};
    PBm = M{3};
    PBp = M{4};
    PC = M{5};
    PD = M{6};
    Pr = M{7};
    PG = M{8};
    
    % Discretization Defect
    delta(:,k) = Reference.x(:,k+1)-xf;
    
    % Compute Discrete Time Dynamics Matrices
    DiscreteDynamics.A(:,:,k) = Pphi;
    DiscreteDynamics.Bm(:,:,k) = Pphi*PBm;
    DiscreteDynamics.Bp(:,:,k) = Pphi*PBp;
    DiscreteDynamics.C(:,:,k) = Pphi*PC;
    DiscreteDynamics.D(:,:,k) = Pphi*PD;
    DiscreteDynamics.r(:,k) = Pphi*Pr;
    DiscreteDynamics.G(:,:,k) = Pphi*PG*Pphi';

    % Append Trajectory:
    DiscreteDynamics.Trajectory = [DiscreteDynamics.Trajectory Trajectory];

end

DiscreteDynamics.E = LinearizedDynamics.E();

end

function [Pdot,Trajectory] = Derivatives(SCP,P,t,up,um,p,h,c,tp,tm,shape,LinearizedDynamics)

% Unvectorize Integration Vector
M = Unflatten(P,shape);
x = M{1};
Pphi = M{2};

% Compute Time Ratios in Integration for First Order Hold
if isequal(SCP.Discretization,'ZOH')
    lambdakm = 1;
    lambdakp = 0;
elseif isequal(SCP.Discretization,'FOH')
    lambdakm = (tp-t)/(tp-tm);
    lambdakp = 1-lambdakm;
else
    error('Invalid Discretization Method Specified (FOH or ZOH)');
end

% Determine Input from First Order Hold
u = lambdakm*um + lambdakp*up;

% Compute Dynamics
f = LinearizedDynamics.f(x,u,p,h,c);
A = LinearizedDynamics.A(x,u,p,h,c);
B = LinearizedDynamics.B(x,u,p,h,c);
C = LinearizedDynamics.C(x,u,p,h,c);
D = LinearizedDynamics.D(x,u,p,h,c);
r = f - A*x - B*u - C*p - D*h; % Residual
G = LinearizedDynamics.g(x,u,p,h,c);

% Invert State Transition Matrix
psi = Pphi^-1;

% Compute Derivatives 
Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PCdot = psi*C;
PDdot = psi*D;
PRdot = psi*r;
PGdot = psi*G*G'*psi';

% Vectorize Result
Pdot = Flatten({Pxdot,Pphidot,PBmdot,PBpdot,PCdot,PDdot,PRdot,PGdot},shape);

% Trajecory Structure:
Trajectory.x = x;
Trajectory.u = u;

end

function [A] = Unflatten(V,shape)

if length(V) ~= sum(prod(shape,2))
    error('Flattening Array Dimensions Incommensurate');
end

ind = 1;
for i = 1:size(shape,1)

    A{i} = reshape(V(ind:ind+prod(shape(i,:))-1),shape(i,:));
    ind = ind+prod(shape(i,:));

end

end

function V = Flatten(A,shape)

if length(A) ~= size(shape,1)
    error('Flattening Array Dimensions Incommensurate');
end

V = [];

for i = 1:size(shape,1)

    if ~all(size(A{i}) == shape(i,:))
        error('Flattening Array Dimensions Incommensurate');
    end

    V = [V;reshape(A{i},[prod(shape(i,:)),1])];

end

end

function [x,VarsOut] = RK4(xdot,tspan,x0)

% Set Initial Condition
x = x0;

% Time Step
dt = tspan(2)-tspan(1);

% RK4
for i = 1:length(tspan)-1
    
    [f1,VarsOut(i)] = xdot(tspan(i),x);
    [f2,~] = xdot(tspan(i)+dt/2,x+f1*dt/2);
    [f3,~] = xdot(tspan(i)+dt/2,x+f2*dt/2);
    [f4,~] = xdot(tspan(i)+dt,x+f3*dt);

    x = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

end

function [DiscreteConstraints] = DiscretizeNonConvexConstraints(Reference,LinearizedConstraints,dims)

% Preallocate
DiscreteConstraints.s = zeros(dims.ns,dims.N);
DiscreteConstraints.As = zeros(dims.ns,dims.nx,dims.N);
DiscreteConstraints.Bs = zeros(dims.ns,dims.nu,dims.N);
DiscreteConstraints.Cs = zeros(dims.ns,dims.np,dims.N);
DiscreteConstraints.Ds = zeros(dims.ns,dims.nh,dims.N);
DiscreteConstraints.rs = zeros(dims.ns,dims.N);

% Discretize
for k = 1:dims.N

    x = Reference.x(:,k);
    u = Reference.u(:,k);
    p = Reference.p;
    h = Reference.h(:,k);
    c = Reference.c;

    DiscreteConstraints.s(:,k) = LinearizedConstraints.s(x,u,p,h,c);
    DiscreteConstraints.As(:,:,k) = LinearizedConstraints.As(x,u,p,h,c);
    DiscreteConstraints.Bs(:,:,k) = LinearizedConstraints.Bs(x,u,p,h,c);
    DiscreteConstraints.Cs(:,:,k) = LinearizedConstraints.Cs(x,u,p,h,c);
    DiscreteConstraints.Ds(:,:,k) = LinearizedConstraints.Ds(x,u,p,h,c);
    DiscreteConstraints.rs(:,k) = DiscreteConstraints.s(:,k) - DiscreteConstraints.As(:,:,k)*x - DiscreteConstraints.Bs(:,:,k)*u - DiscreteConstraints.Cs(:,:,k)*p - DiscreteConstraints.Ds(:,:,k)*h;

end

end

function [DiscreteCost] = DiscretizeCost(Reference,LinearizedCost,dims)

% Terminal Cost:
DiscreteCost.phi = LinearizedCost.phi(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteCost.Aphi = LinearizedCost.Aphi(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteCost.Bphi = LinearizedCost.Bphi(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteCost.Cphi = LinearizedCost.Cphi(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteCost.Dphi = LinearizedCost.Dphi(Reference.x(:,end),Reference.u(:,end),Reference.p,Reference.h(:,end),Reference.c);
DiscreteCost.rphi = DiscreteCost.phi - DiscreteCost.Aphi*Reference.x(:,end)-DiscreteCost.Bphi*Reference.u(:,end)-DiscreteCost.Cphi*Reference.p-DiscreteCost.Dphi*Reference.h(:,end);

% Preallocate
DiscreteCost.L = zeros(1,dims.N);
DiscreteCost.AL = zeros(1,dims.nx,dims.N);
DiscreteCost.BL = zeros(1,dims.nu,dims.N);
DiscreteCost.CL = zeros(1,dims.np,dims.N);
DiscreteCost.DL = zeros(1,dims.nh,dims.N);
DiscreteCost.rL = zeros(1,dims.N);

% Discretize
for k = 1:dims.N

    x = Reference.x(:,k);
    u = Reference.u(:,k);
    p = Reference.p;
    h = Reference.h(:,k);
    c = Reference.c;

    DiscreteCost.L(:,k) = LinearizedCost.L(x,u,p,h,c);
    DiscreteCost.AL(:,:,k) = LinearizedCost.AL(x,u,p,h,c);
    DiscreteCost.BL(:,:,k) = LinearizedCost.BL(x,u,p,h,c);
    DiscreteCost.CL(:,:,k) = LinearizedCost.CL(x,u,p,h,c);
    DiscreteCost.DL(:,:,k) = LinearizedCost.DL(x,u,p,h,c);
    DiscreteCost.rL(:,k) = DiscreteCost.L(:,k) - DiscreteCost.AL(:,:,k)*x - DiscreteCost.BL(:,:,k)*u - DiscreteCost.CL(:,:,k)*p - DiscreteCost.DL(:,:,k)*h;

end

end