function [LinearizedCost,LinearizedDynamics,LinearizedConstraints] = Linearize(SCP,System,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition)

syms t; % Time Variable
[x,u,c] = Variables(System);
[f,slack] = Dynamics(System,x,u,c);
[s] = NonConvexConstraints(System,x,u,c);
[g0] = InitialCondition(System,x,c);
[gf] = FinalCondition(System,x,c);
[phi,L] = Cost(System,t,x,u,c);

% Dynamics Jacobians
if ~SCP.FreeTime && ~SCP.AdaptiveMesh
    LinearizedDynamics.f = matlabFunction(System.T*f,"Vars",{x,u,c});
    LinearizedDynamics.A = matlabFunction(System.T*Jacobian(f,x),"Vars",{x,u,c});
    LinearizedDynamics.B = matlabFunction(System.T*Jacobian(f,u),"Vars",{x,u,c});
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

if ~isempty(slack)
    LinearizedDynamics.E = matlabFunction(Jacobian(x,slack));
else
    LinearizedDynamics.E = zeros(length(x),0);
end

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

% Cost Jacobians
LinearizedCost.phi = matlabFunction(phi,"Vars",{t,x,c});
LinearizedCost.Aphi = matlabFunction(Jacobian(phi,x),"Vars",{t,x,c});
LinearizedCost.Fphi = matlabFunction(Jacobian(phi,t),"Vars",{t,x,c});

LinearizedCost.L = matlabFunction(L,"Vars",{t,x,u,c});
LinearizedCost.AL = matlabFunction(Jacobian(L,x),"Vars",{t,x,u,c});
LinearizedCost.BL = matlabFunction(Jacobian(L,u),"Vars",{t,x,u,c});
LinearizedCost.FL = matlabFunction(Jacobian(L,t),"Vars",{t,x,u,c});

end

function fx = Jacobian(f,x)

if isempty(f)
    fx = sym(zeros(0,length(x)));
    return;
end

fx = [];

for i = 1:length(x)
    fx = [fx diff(f,x(i))];
end

end