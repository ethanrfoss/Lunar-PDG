function [DiscreteCost,DiscreteDynamics,DiscreteConstraints,delta] = Discretize(SCP,System,Reference,LinearizedCost,LinearizedDynamics,LinearizedConstraints,dims)

if SCP.FreeTime && ~SCP.AdaptiveMesh
    [DiscreteDynamics,delta] = DiscretizeDynamicsFreeTime(SCP,Reference,LinearizedDynamics,dims);
    T = Reference.p;
end
if ~SCP.FreeTime && ~SCP.AdaptiveMesh
    [DiscreteDynamics,delta] = DiscretizeDynamicsFixedTime(SCP,Reference,LinearizedDynamics,dims);
    T = System.T;
end
if SCP.AdaptiveMesh
    [DiscreteDynamics,delta] = DiscretizeDynamicsAdaptiveTime(SCP,Reference,LinearizedDynamics,dims);
    T = sum(Reference.h);
end

if dims.ns ~= 0
    DiscreteConstraints = DiscretizeNonConvexConstraints(Reference,LinearizedConstraints,dims);
end

DiscreteConstraints.H0 = LinearizedConstraints.H0(Reference.x(:,1),Reference.c);
DiscreteConstraints.l0 = LinearizedConstraints.g0(Reference.x(:,1),Reference.c) - DiscreteConstraints.H0*Reference.x(:,1);
DiscreteConstraints.Hf = LinearizedConstraints.Hf(Reference.x(:,end),Reference.c);
DiscreteConstraints.lf = LinearizedConstraints.gf(Reference.x(:,end),Reference.c) - DiscreteConstraints.Hf*Reference.x(:,end);

[DiscreteCost] = DiscretizeCost(Reference,T,LinearizedCost,dims);

end

function [DiscreteDynamics,delta] = DiscretizeDynamicsFreeTime(SCP,Reference,LinearizedDynamics,dims)

% Preallocate:
delta = zeros(dims.nx,dims.N-1);
DiscreteDynamics.A = zeros(dims.nx,dims.nx,dims.N-1);
DiscreteDynamics.Bm = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.Bp = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.F = zeros(dims.nx,1,dims.N-1);
DiscreteDynamics.r = zeros(dims.nx,dims.N-1);

% Flattening Shape:
shape = [dims.nx,1;
         dims.nx,dims.nx;
         dims.nx,dims.nu;
         dims.nx,dims.nu;
         dims.nx,1;
         dims.nx,1];

% Time:
t = linspace(0,1,dims.N);

% For Loop to Compute Discrete Time Dynamics at Each Time Step
for k = 1:dims.N-1
    
    % Initial Conditions
    phi0 = eye(dims.nx); % Initial Condition for A
    PBm0 = zeros(dims.nx,dims.nu); % Initial Condition for Bm
    PBp0 = zeros(dims.nx,dims.nu); % Initial Condition for Bp
    PF0 = zeros(dims.nx,1);  % Initial Condition for F
    Pr0 = zeros(dims.nx,1); % Initial Condition for r

    % Vectorize Initial Conditions
    P0 = Flatten({Reference.x(:,k),phi0,PBm0,PBp0,PF0,Pr0},shape);
    
    % Integrate from current time step to next:
    tspan = linspace(t(k),t(k+1),SCP.Nsub);
    P = RK4(@(tc,x) DerivativesFreeTime(x,tc,Reference.u(:,k+1),Reference.u(:,k),Reference.p,Reference.c,t(k+1),t(k),shape,LinearizedDynamics),tspan,P0);
    
    % Unvectorize Integrated Result
    M = Unflatten(P,shape);
    xf = M{1};
    Pphi = M{2};
    PBm = M{3};
    PBp = M{4};
    PF = M{5};
    Pr = M{6};
    
    % Discretization Defect
    delta(:,k) = Reference.x(:,k+1)-xf;
    
    % Compute Discrete Time Dynamics Matrices
    DiscreteDynamics.A(:,:,k) = Pphi;
    DiscreteDynamics.Bm(:,:,k) = Pphi*PBm;
    DiscreteDynamics.Bp(:,:,k) = Pphi*PBp;
    DiscreteDynamics.F(:,:,k) = Pphi*PF;
    DiscreteDynamics.r(:,k) = Pphi*Pr;

end

DiscreteDynamics.E = LinearizedDynamics.E();

end

function [Pdot] = DerivativesFreeTime(P,t,up,um,p,c,tp,tm,shape,LinearizedDynamics)

% Unvectorize Integration Vector
M = Unflatten(P,shape);
x = M{1};
Pphi = M{2};

% Compute Time Ratios in Integration for First Order Hold
lambdakm = (tp-t)/(tp-tm);
lambdakp = 1-lambdakm;

% Determine Input from First Order Hold
u = lambdakm*um + lambdakp*up;

% Compute Dynamics
f = LinearizedDynamics.f(x,u,p,c);
A = LinearizedDynamics.A(x,u,p,c);
B = LinearizedDynamics.B(x,u,p,c);
F = LinearizedDynamics.F(x,u,p,c);
r = f - A*x - B*u - F*p; % Residual

% Invert State Transition Matrix
psi = Pphi^-1;

% Compute Derivatives 
Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PFdot = psi*F;
PRdot = psi*r;

% Vectorize Result
Pdot = Flatten({Pxdot,Pphidot,PBmdot,PBpdot,PFdot,PRdot},shape);

end

function [DiscreteDynamics,delta] = DiscretizeDynamicsFixedTime(SCP,Reference,LinearizedDynamics,dims)

% Preallocate:
delta = zeros(dims.nx,dims.N-1);
DiscreteDynamics.A = zeros(dims.nx,dims.nx,dims.N-1);
DiscreteDynamics.Bm = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.Bp = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.r = zeros(dims.nx,dims.N-1);

% Flattening Shape:
shape = [dims.nx,1;
         dims.nx,dims.nx;
         dims.nx,dims.nu;
         dims.nx,dims.nu;
         dims.nx,1];

% Time:
t = linspace(0,1,dims.N);

% For Loop to Compute Discrete Time Dynamics at Each Time Step
for k = 1:dims.N-1
    
    % Initial Conditions
    phi0 = eye(dims.nx); % Initial Condition for A
    PBm0 = zeros(dims.nx,dims.nu); % Initial Condition for Bm
    PBp0 = zeros(dims.nx,dims.nu); % Initial Condition for Bp
    Pr0 = zeros(dims.nx,1); % Initial Condition for r

    % Vectorize Initial Conditions
    P0 = Flatten({Reference.x(:,k),phi0,PBm0,PBp0,Pr0},shape);
    
    % Integrate from current time step to next:
    tspan = linspace(t(k),t(k+1),SCP.Nsub);
    P = RK4(@(tc,x) DerivativesFixedTime(x,tc,Reference.u(:,k+1),Reference.u(:,k),Reference.c,t(k+1),t(k),shape,LinearizedDynamics),tspan,P0);
    
    % Unvectorize Integrated Result
    M = Unflatten(P,shape);
    xf = M{1};
    Pphi = M{2};
    PBm = M{3};
    PBp = M{4};
    Pr = M{5};
    
    % Discretization Defect
    delta(1:dims.nx,k) = Reference.x(:,k+1)-xf;
    
    % Compute Discrete Time Dynamics Matrices
    DiscreteDynamics.A(:,:,k) = Pphi;
    DiscreteDynamics.Bm(:,:,k) = Pphi*PBm;
    DiscreteDynamics.Bp(:,:,k) = Pphi*PBp;
    DiscreteDynamics.r(:,k) = Pphi*Pr;

end

DiscreteDynamics.E = LinearizedDynamics.E();

end

function [Pdot] = DerivativesFixedTime(P,t,up,um,c,tp,tm,shape,LinearizedDynamics)

% Unvectorize Integration Vector
M = Unflatten(P,shape);
x = M{1};
Pphi = M{2};

% Compute Time Ratios in Integration for First Order Hold
lambdakm = (tp-t)/(tp-tm);
lambdakp = 1-lambdakm;

% Determine Input from First Order Hold
u = lambdakm*um + lambdakp*up;

% Compute Dynamics
f = LinearizedDynamics.f(x,u,c);
A = LinearizedDynamics.A(x,u,c);
B = LinearizedDynamics.B(x,u,c);
r = f - A*x - B*u; % Residual

% Invert State Transition Matrix
psi = Pphi^-1;

% Compute Derivatives 
Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PRdot = psi*r;

% Vectorize Result
Pdot = Flatten({Pxdot,Pphidot,PBmdot,PBpdot,PRdot},shape);

end

function [DiscreteDynamics,delta] = DiscretizeDynamicsAdaptiveTime(SCP,Reference,LinearizedDynamics,dims)

% Preallocate:
delta = zeros(dims.nx,dims.N-1);
DiscreteDynamics.A = zeros(dims.nx,dims.nx,dims.N-1);
DiscreteDynamics.Bm = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.Bp = zeros(dims.nx,dims.nu,dims.N-1);
DiscreteDynamics.H = zeros(dims.nx,1,dims.N-1);
DiscreteDynamics.r = zeros(dims.nx,dims.N-1);

% Flattening Shape:
shape = [dims.nx,1;
         dims.nx,dims.nx;
         dims.nx,dims.nu;
         dims.nx,dims.nu;
         dims.nx,1;
         dims.nx,1];

% Time:
t = linspace(0,1,dims.N);

% For Loop to Compute Discrete Time Dynamics at Each Time Step
for k = 1:dims.N-1
    
    % Initial Conditions
    phi0 = eye(dims.nx); % Initial Condition for A
    PBm0 = zeros(dims.nx,dims.nu); % Initial Condition for Bm
    PBp0 = zeros(dims.nx,dims.nu); % Initial Condition for Bp
    PH0 = zeros(dims.nx,1);  % Initial Condition for F
    Pr0 = zeros(dims.nx,1); % Initial Condition for r

    % Vectorize Initial Conditions
    P0 = Flatten({Reference.x(:,k),phi0,PBm0,PBp0,PH0,Pr0},shape);
    
    % Integrate from current time step to next:
    tspan = linspace(t(k),t(k+1),SCP.Nsub);
    P = RK4(@(tc,x) DerivativesAdaptiveTime(x,tc,Reference.u(:,k+1),Reference.u(:,k),Reference.h(k),Reference.c,t(k+1),t(k),shape,LinearizedDynamics),tspan,P0);
    
    % Unvectorize Integrated Result
    M = Unflatten(P,shape);
    xf = M{1};
    Pphi = M{2};
    PBm = M{3};
    PBp = M{4};
    PH = M{5};
    Pr = M{6};
    
    % Discretization Defect
    delta(1:dims.nx,k) = Reference.x(:,k+1)-xf;
    
    % Compute Discrete Time Dynamics Matrices
    DiscreteDynamics.A(:,:,k) = Pphi;
    DiscreteDynamics.Bm(:,:,k) = Pphi*PBm;
    DiscreteDynamics.Bp(:,:,k) = Pphi*PBp;
    DiscreteDynamics.H(:,:,k) = Pphi*PH;
    DiscreteDynamics.r(:,k) = Pphi*Pr;

end

DiscreteDynamics.E = LinearizedDynamics.E();

end

function [Pdot] = DerivativesAdaptiveTime(P,t,up,um,h,c,tp,tm,shape,LinearizedDynamics)

% Unvectorize Integration Vector
M = Unflatten(P,shape);
x = M{1};
Pphi = M{2};

% Compute Time Ratios in Integration for First Order Hold
lambdakm = (tp-t)/(tp-tm);
lambdakp = 1-lambdakm;

% Determine Input from First Order Hold
u = lambdakm*um + lambdakp*up;

% Compute Dynamics
f = LinearizedDynamics.f(x,u,h,c);
A = LinearizedDynamics.A(x,u,h,c);
B = LinearizedDynamics.B(x,u,h,c);
H = LinearizedDynamics.H(x,u,h,c);
r = f - A*x - B*u - H*h; % Residual

% Invert State Transition Matrix
psi = Pphi^-1;

% Compute Derivatives 
Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PHdot = psi*H;
PRdot = psi*r;

% Vectorize Result
Pdot = Flatten({Pxdot,Pphidot,PBmdot,PBpdot,PHdot,PRdot},shape);

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

function x = RK4(xdot,tspan,x0)

% Set Initial Condition
x = x0;

% Time Step
dt = tspan(2)-tspan(1);

% RK4
for i = 1:length(tspan)-1
    
    f1 = xdot(tspan(i),x);
    f2 = xdot(tspan(i)+dt/2,x+f1*dt/2);
    f3 = xdot(tspan(i)+dt/2,x+f2*dt/2);
    f4 = xdot(tspan(i)+dt,x+f3*dt);

    x = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

end

function [DiscreteConstraints] = DiscretizeNonConvexConstraints(Reference,LinearizedConstraints,dims)

% Preallocate
DiscreteConstraints.s = zeros(dims.ns,dims.N);
DiscreteConstraints.C = zeros(dims.ns,dims.nx,dims.N);
DiscreteConstraints.D = zeros(dims.ns,dims.nu,dims.N);
DiscreteConstraints.rp = zeros(dims.ns,dims.N);

% Discretize
for k = 1:dims.N

    x = Reference.x(:,k);
    u = Reference.u(:,k);
    c = Reference.c;

    DiscreteConstraints.s(:,k) = LinearizedConstraints.s(x,u,c);
    DiscreteConstraints.C(:,:,k) = LinearizedConstraints.C(x,u,c);
    DiscreteConstraints.D(:,:,k) = LinearizedConstraints.D(x,u,c);
    DiscreteConstraints.rp(:,k) = DiscreteConstraints.s(:,k) - DiscreteConstraints.C(:,:,k)*x - DiscreteConstraints.D(:,:,k)*u;

end

end

function [DiscreteCost] = DiscretizeCost(Reference,T,LinearizedCost,dims)

% Terminal Cost:
DiscreteCost.phi = LinearizedCost.phi(T,Reference.x(:,end),Reference.c);
DiscreteCost.Aphi = LinearizedCost.Aphi(T,Reference.x(:,end),Reference.c);
DiscreteCost.Fphi = LinearizedCost.Fphi(T,Reference.x(:,end),Reference.c);
DiscreteCost.rphi = DiscreteCost.phi - DiscreteCost.Aphi*Reference.x(:,end)-DiscreteCost.Fphi*T;

% Preallocate
DiscreteCost.L = zeros(1,dims.N);
DiscreteCost.AL = zeros(1,dims.nx,dims.N);
DiscreteCost.BL = zeros(1,dims.nu,dims.N);
DiscreteCost.FL = zeros(1,1,dims.N);
DiscreteCost.rL = zeros(1,dims.N);

% Discretize
for k = 1:dims.N

    x = Reference.x(:,k);
    u = Reference.u(:,k);
    c = Reference.c;

    DiscreteCost.L(:,k) = LinearizedCost.L(T,x,u,c);
    DiscreteCost.AL(:,:,k) = LinearizedCost.AL(T,x,u,c);
    DiscreteCost.BL(:,:,k) = LinearizedCost.BL(T,x,u,c);
    DiscreteCost.FL(:,:,k) = LinearizedCost.FL(T,x,u,c);
    DiscreteCost.rL(:,k) = DiscreteCost.L(:,k) - DiscreteCost.AL(:,:,k)*x - DiscreteCost.BL(:,:,k)*u - DiscreteCost.FL(:,:,k)*T;

end

end