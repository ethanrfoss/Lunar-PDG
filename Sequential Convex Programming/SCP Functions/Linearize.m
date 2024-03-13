function [LinearizedCost,LinearizedDynamics,LinearizedConstraints] = Linearize(System,Variables,Cost,Dynamics,NonConvexConstraints,InitialCondition,FinalCondition)

[x,u,p,h,c] = Variables(System);
[f,g,slack] = Dynamics(System,x,u,p,h,c);
[s] = NonConvexConstraints(System,x,u,p,h,c);
[g0] = InitialCondition(System,x,u,p,h,c);
[gf] = FinalCondition(System,x,u,p,h,c);
[phi,L] = Cost(System,x,u,p,h,c);

% Make Sure all functions and variables are sym type:
x = sym(x);
u = sym(u);
p = sym(p);
h = sym(h);
c = sym(c);
f = sym(f);
g = sym(g);
slack = sym(slack);
s = sym(s);
g0 = sym(g0);
gf = sym(gf);
phi = sym(phi);
L = sym(L);

% Dynamics Jacobians
LinearizedDynamics.f = matlabFunction(f,"Vars",{x,u,p,h,c});
LinearizedDynamics.A = matlabFunction(Jacobian(f,x),"Vars",{x,u,p,h,c});
LinearizedDynamics.B = matlabFunction(Jacobian(f,u),"Vars",{x,u,p,h,c});
LinearizedDynamics.C = matlabFunction(Jacobian(f,p),"Vars",{x,u,p,h,c});
LinearizedDynamics.D = matlabFunction(Jacobian(f,h),"Vars",{x,u,p,h,c});
LinearizedDynamics.E = matlabFunction(Jacobian(x,slack));
LinearizedDynamics.g = matlabFunction(g,"Vars",{x,u,p,h,c});

% Non Convex Constraints Jacobians
LinearizedConstraints.s = matlabFunction(s,"Vars",{x,u,p,h,c});
LinearizedConstraints.As = matlabFunction(Jacobian(s,x),"Vars",{x,u,p,h,c});
LinearizedConstraints.Bs = matlabFunction(Jacobian(s,u),"Vars",{x,u,p,h,c});
LinearizedConstraints.Cs = matlabFunction(Jacobian(s,p),"Vars",{x,u,p,h,c});
LinearizedConstraints.Ds = matlabFunction(Jacobian(s,h),"Vars",{x,u,p,h,c});

% Initial Condition Jacobians
LinearizedConstraints.g0 = matlabFunction(g0,"Vars",{x,u,p,h,c});
LinearizedConstraints.H0 = matlabFunction(Jacobian(g0,x),"Vars",{x,u,p,h,c});
LinearizedConstraints.I0 = matlabFunction(Jacobian(g0,u),"Vars",{x,u,p,h,c});
LinearizedConstraints.J0 = matlabFunction(Jacobian(g0,p),"Vars",{x,u,p,h,c});
LinearizedConstraints.K0 = matlabFunction(Jacobian(g0,h),"Vars",{x,u,p,h,c});

% Final Condition Jacobians
LinearizedConstraints.gf = matlabFunction(gf,"Vars",{x,u,p,h,c});
LinearizedConstraints.Hf = matlabFunction(Jacobian(gf,x),"Vars",{x,u,p,h,c});
LinearizedConstraints.If = matlabFunction(Jacobian(gf,u),"Vars",{x,u,p,h,c});
LinearizedConstraints.Jf = matlabFunction(Jacobian(gf,p),"Vars",{x,u,p,h,c});
LinearizedConstraints.Kf = matlabFunction(Jacobian(gf,h),"Vars",{x,u,p,h,c});

% Cost Jacobians
LinearizedCost.phi = matlabFunction(phi,"Vars",{x,u,p,h,c});
LinearizedCost.Aphi = matlabFunction(Jacobian(phi,x),"Vars",{x,u,p,h,c});
LinearizedCost.Bphi = matlabFunction(Jacobian(phi,u),"Vars",{x,u,p,h,c});
LinearizedCost.Cphi = matlabFunction(Jacobian(phi,p),"Vars",{x,u,p,h,c});
LinearizedCost.Dphi = matlabFunction(Jacobian(phi,h),"Vars",{x,u,p,h,c});

LinearizedCost.L = matlabFunction(L,"Vars",{x,u,p,h,c});
LinearizedCost.AL = matlabFunction(Jacobian(L,x),"Vars",{x,u,p,h,c});
LinearizedCost.BL = matlabFunction(Jacobian(L,u),"Vars",{x,u,p,h,c});
LinearizedCost.CL = matlabFunction(Jacobian(L,p),"Vars",{x,u,p,h,c});
LinearizedCost.DL = matlabFunction(Jacobian(L,h),"Vars",{x,u,p,h,c});

end

function fx = Jacobian(f,x)

if isempty(f)
    fx = sym(zeros(0,length(x)));
    return;
end

if isempty(x)
    fx = sym(zeros(length(f),0));
    return;
end

fx = [];

for i = 1:length(x)
    fx = [fx diff(f,x(i))];
end

end