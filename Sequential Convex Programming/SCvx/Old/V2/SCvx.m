
function [L,x,u,p,h,vars] = SCvx(SCP,System,xbar,ubar,pbar,hbar,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition,Scaling,Continuation)

% Add Solver Path:
if strcmp(SCP.Solver,'ECOS')
    addpath([fileparts(mfilename('fullpath')) '\ecos-matlab-master']);
end
if strcmp(SCP.Solver,'MOSEK')
    addpath('C:\Program Files\Mosek\10.1\toolbox\r2017a');
end

% Algorithm Start
tic;

% Check Dimension Commensurabilites and Get Dimensions:
[N,nx,nu,np,nh,nv,ns,nic,ntc,nlx,nlu,nlp,nlh,nlht,nqx,nqu,nqp,nbx,nbu] = Dimensions(System,xbar,ubar,pbar,hbar,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition);

% Time Vector:
t = linspace(0,1,N);

% Cost:
[lx,lu,lp,Bx,Bu,Bh] = Cost();

% Scaling Matrices:
[Sx,Su,Sp,Sh,cx,cu,cp,ch] = ScalingMatrices(System,Scaling);

% Discretize:
[A,Bm,Bp,F,H,r,E,delta] = DiscretizeDynamics(SCP,System,xbar,ubar,pbar,hbar,t,nx,nu,np,nv,N,Dynamics);
[WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp] = DiscretizeConvexConstraint(System,xbar,ubar,pbar,hbar,t,nx,nu,nh,nlx,nlu,nlh,nlht,nqx,nqu,N,ConvexConstraints);
[s,C,D,G,I,rp] = DiscretizeNonConvexConstraints(System,xbar,ubar,pbar,hbar,t,nx,nu,np,nh,ns,N,NonConvexConstraints);
[gic,H0,K0] = InitialCondition(System,xbar(:,1),pbar); l0 = gic - H0*xbar(:,1) - K0*pbar;
[gtc,Hf,Kf] = FinalCondition(System,xbar(:,end),pbar); lf = gtc - Hf*xbar(:,end) - Kf*pbar;

% Initial NonLiear Cost
Jbar = NonLinearCost(SCP,xbar,ubar,pbar,hbar,delta,N,lx,lu,lp,Bx,Bu,Bh,s,gic,gtc);

% Initialize Continuation Parameter:
gamma = 0;

iter = 1;

while true
    
    % Solve Problem as a SOCP
    if strcmp(SCP.Solver,'ECOS')
        [cs,As,bs,Gs,hs,Ls,Qs] = Parse(SCP,xbar,ubar,pbar,hbar,lx,lu,lp,Bx,Bu,Bh,A,Bm,Bp,F,H,r,E,WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp,C,D,G,rp,H0,K0,l0,Hf,Kf,lf,Sx,Su,Sp,cx,cu,cp,N,nx,nu,np,nh,nv,ns,nic,ntc,nlx,nlu,nlp,nlh,nlht,nqx,nqu,nqp,nbx,nbu);
        [Lhat,z] = ConvexSolverECOS(cs,As,bs,Gs,hs,Ls,Qs);
        [Lstar,xstar,ustar,pstar,vars(iter)] = Descale(Lhat,z,Sx,cx,Su,cu,Sp,cp,lx,lu,lp,N,nx,nu,np,nv,ns,nic,ntc);
    elseif strcmp(SCP.Solver,'MOSEK')
        [cs,As,bs,Gs,hs,Ls,Qs] = Parse(SCP,xbar,ubar,pbar,lx,lu,lp,Bx,Bu,A,Bm,Bp,F,r,E,WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp,C,D,G,rp,H0,K0,l0,Hf,Kf,lf,Sx,Su,Sp,cx,cu,cp,N,nx,nu,np,nv,ns,nic,ntc,nlx,nlu,nlp,nqx,nqu,nqp,nbx,nbu);
        [Lhat,z] = ConvexSolverMOSEK(cs,As,bs,Gs,hs,Ls,Qs);
        [Lstar,xstar,ustar,pstar,vars(iter)] = Descale(Lhat,z,Sx,cx,Su,cu,Sp,cp,lx,lu,lp,N,nx,nu,np,nv,ns,nic,ntc);
    elseif strcmp(SCP.Solver,'CVX')
        [Lstar,xstar,ustar,pstar,vars(iter)] = ConvexSolverCVX(SCP,xbar,ubar,pbar,lx,lu,lp,Bx,Bu,A,Bm,Bp,F,r,E,WLx,yLx,WLu,yLu,WLp,yLp,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp,C,D,G,rp,H0,K0,l0,Hf,Kf,lf,Sx,Su,Sp,cx,cu,cp,N,nx,nu,np,nv,ns,nic,ntc,nlx,nlu,nlp);
    else
        disp('Invalid Convex Solver Specified');
    end

    % Check Convergence Criterion
    if abs(Jbar-Lstar) <= SCP.epsilon && gamma == 1 %norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)) <= SCP.epsilon
        disp(['SCP Algorithm Converged in ' num2str(toc) ' seconds']);
        break;
    else
        disp(['Current Iteration: ' num2str(iter) ' Current Tolerance: ' num2str(abs(Jbar-Lstar))]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
    end

    if iter >= SCP.maxIter
        disp(['Maximum Iteration Count Reached. Current Tolerance: ' num2str(abs(Jbar-Lstar))]); %num2str(norm(pstar-pbar,2) + max(vecnorm(xstar-xbar)))]);
        break
    end
    
    % Redescretize to determine new non linear cost
    [Astar,Bmstar,Bpstar,Fstar,rstar,Estar,delta] = DiscretizeDynamics(SCP,System,xstar,ustar,pstar,t,nx,nu,np,nv,N,Dynamics);
    [sstar,Cstar,Dstar,Gstar,rpstar] = DiscretizeNonConvexConstraints(System,xstar,ustar,pstar,t,nx,nu,np,ns,N,NonConvexConstraints);
    [gicstar,H0star,K0star] = InitialCondition(System,xstar(:,1),pstar); l0star = gicstar - H0star*xstar(:,1) - K0star*pstar;
    [gtcstar,Hfstar,Kfstar] = FinalCondition(System,xstar(:,end),pstar); lfstar = gtcstar - Hfstar*xstar(:,end) - Kfstar*pstar;

    % New Nonlinear Cost
    Jstar = NonLinearCost(SCP,xstar,ustar,pstar,delta,N,lx,lu,lp,Bx,Bu,sstar,gicstar,gtcstar);
    
    % Calculate Convexification Accuracy
    rho = (Jbar - Jstar)/(Jbar-Lstar);
    
    % Update trust region
     if rho < SCP.rho0
        % Update eta, dont update trajectory
        SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
        SCP.etap = max(SCP.etap0,SCP.etap/SCP.betash);
    else
        if rho < SCP.rho1
            % Update eta
            SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
            SCP.etap = max(SCP.etap0,SCP.etap/SCP.betash);
        end
        % If rho1 < rho <= rho2, do not update eta
        if rho >= SCP.rho2
            % Update eta
            SCP.eta = min(SCP.eta1,SCP.betagr*SCP.eta);
            SCP.etap = min(SCP.etap1,SCP.betagr*SCP.etap);
        end
        % Update Nominal Solution
        xbar = xstar;
        ubar = ustar;
        pbar = pstar;
        % Update Nonlinear Cost
        Jbar = Jstar;
        % Update Dynamics
        A = Astar; Bm = Bmstar; Bp = Bpstar;F = Fstar; E = Estar; r = rstar;
        % Update NonConvex Constraint:
        C = Cstar; D = Dstar; G = Gstar; rp = rpstar;
        % Update Initial Condition:
        gic = gicstar; H0 = H0star; K0 = K0star; l0 = l0star;
        % Update Final Condition:
        gtc = gtcstar; Hf = Hfstar; Kf = Kfstar; lf = lfstar;
        % Discretize Rest of Problem:
        [WLx,yLx,WLu,yLu,WLp,yLp,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp] = DiscretizeConvexConstraint(System,xbar,ubar,pbar,t,nx,nu,nlx,nlu,nqx,nqu,N,ConvexConstraints);
    end
    % Update Result:
    L(iter) = Lstar;
    x(1:nx,1:N,iter) = xstar;
    u(1:nu,1:N,iter) = ustar;
    p(1:np,iter) = pstar;

    % Apply Continuation:
    if nargin > 12
        if max(vecnorm(delta)) <= SCP.gammaeps && gamma < 1
            gamma = min(gamma + SCP.dgamma,1);
            disp(['Apply Continuation. gamma = ' num2str(gamma)]);
            [System] = Continuation(System,gamma);
        end
    end

    iter = iter + 1;
    
end

end

%% Dimensions
% Verifies the commensurability of the dimensions of the matrices specified
% for the SCP. Also returns important dimensions
function [N,nx,nu,np,nh,nv,ns,nic,ntc,nlx,nlu,nlp,nlh,nlht,nqx,nqu,nqp,nbx,nbu] = Dimensions(System,xbar,ubar,pbar,hbar,Cost,Dynamics,ConvexConstraints,NonConvexConstraints,InitialCondition,FinalCondition)

[lx,lu,lp,Bx,Bu,Bh] = Cost();
[f,A,B,F,H,E] = Dynamics(System,xbar(:,1),ubar(:,1),pbar,hbar(:,1),0);
[WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp] = ConvexConstraints(System,xbar(:,1),ubar(:,1),pbar,hbar(:,1),0);
[s,C,D,G,I] = NonConvexConstraints(System,xbar(:,1),ubar(:,1),pbar,hbar(:,1),0);
[gic,H0,K0] = InitialCondition(System,xbar(:,1),pbar);
[gtc,Hf,Kf] = FinalCondition(System,xbar(:,end),pbar);

% Number of Nodes:
N = size(xbar,2);
if size(ubar,2) ~= N  || size(hbar,2) ~= N
    error('Incommensurate Node Counts Specified');
end

% State Dimension:
nx = size(xbar,1); 
% if size(lx,1) ~= nx || size(Bx,2) ~= nx || size(f,1) ~= nx || size(A,1) ~= nx || size(A,2) ~= nx || size(B,1) ~= nx || size(F,1) ~= nx || size(E,1) ~= nx || size(WLx,2) ~= nx  || size(WQx,2) ~= nx || size(C,2) ~= nx || size(H0,2) ~= nx || size(Hf,2) ~= nx
%     error('Incommensurate State Dimensions Specified');
% end

% Input Dimension:
nu = size(ubar,1); 
% if size(lu,1) ~= nu || size(Bu,2) ~= nu || size(B,2) ~= nu || size(WLu,2) ~= nu || size(WQu,2) ~= nu || size(D,2) ~= nu
%     error('Incommensurate Input Dimensions Specified');
% end

% Parameter Dimension:
np = size(pbar,2);
% if size(lp,1) ~= np || size(F,2) ~= np || size(WLp,2) ~= np || size(WQp,2) ~= np || size(G,2) ~= np || size(K0,2) ~= np || size(Kf,2) ~= np
%     error('Incommensurate Parameter Dimensions Specified');
% end

% Time Varying Parameter Dimension:
nh = size(hbar,2);

% Dynamic Slack Dimension:
nv = size(E,2);

% Constraint Slack Dimension:
ns = size(s,1);
% if size(C,2) ~= ns || size(D,2) ~= ns || size(G,2) ~= ns
%     error('Incommensurate Non Convex Constraint Dimensions Specified');
% end

% Initial Condition Dimension:
nic = size(gic,1);
if size(H0,1) ~= nic || size(K0,1) ~= nic
    error('Incommensurate Initial Condition Dimensions Specified');
end

% Final Condition Dimension:
ntc = size(gtc,1);
if size(Hf,1) ~= ntc || size(Kf,1) ~= ntc
    error('Incommensurate Final Condition Dimensions Specified');
end

% Convex Linear Cone State Dimensions
nlx = size(WLx,1);
if size(yLx,1) ~= nlx
    error('Incommensurate Covex Linear Cone State Dimensions Specified');
end

% Convex Linear Cone Input Dimensions
nlu = size(WLu,1);
if size(yLu,1) ~= nlu
    error('Incommensurate Covex Linear Cone Input Dimensions Specified');
end

% Convex Linear Cone Parameter Dimensions
nlp = size(WLp,1);
if size(yLp,1) ~= nlp
    error('Incommensurate Covex Linear Cone Parameter Dimensions Specified');
end

% Convex Linear Cone Time Varying Parameter Dimensions
nlh = size(WLh,1);
if size(yLh,1) ~= nlh
    error('Incommensurate Covex Linear Cone Time VaryingParameter Dimensions Specified');
end

% Convex Linear Cone Time Varying Parameter Dimensions
nlht = size(WLht,1);
if size(yLht,1) ~= nlht
    error('Incommensurate Covex Linear Cone Time Varying Parameter Dimensions Specified');
end

% Convex Quadratic Cone State Dimensions
nqx = size(WQx,1);
if size(yQx,1) ~= nqx || sum(dqx) ~= nqx
    error('Incommensurate Covex Quadratic Cone State Dimensions Specified');
end

% Convex Quadratic Cone Input Dimensions
nqu = size(WQu,1);
if size(yQu,1) ~= nqu || sum(dqu) ~= nqu
    error('Incommensurate Covex Quadratic Cone Input Dimensions Specified');
end

% Convex Quadratic Cone Parameter Dimensions
nqp = size(WQp,1);
if size(yQp,1) ~= nqp || sum(dqp) ~= nqp
    error('Incommensurate Covex Quadratic Cone Parameter Dimensions Specified');
end

% Integral Induced Norm State Cost Dimension
nbx = size(Bx,1);

% Integral Induced Norm Input Cost Dimension
nbu = size(Bu,1);

end

function [Lhat,z] = ConvexSolverECOS(c,A,b,G,h,L,Q)

% Set Linear and Quadratic Cone Dimensions
dims.l = L;
dims.q = Q;

% Call Solver
[z,~,info] = ecos(c,sparse(G),h,dims,sparse(A),b);

Lhat  = info.pcost; % Scaled Cost

end

function [L,x,u,p,vars] = ConvexSolverCVX(SCP,xbar,ubar,pbar,lx,lu,lp,Bx,Bu,A,Bm,Bp,F,r,E,WLx,yLx,WLu,yLu,WLp,yLp,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp,C,D,G,rp,H0,K0,l0,Hf,Kf,lf,Sx,Su,Sp,cx,cu,cp,N,nx,nu,np,nv,ns,nic,ntc,nlx,nlu,nlp)

cvx_begin

variables x(nx,N) u(nu,N) p(np) v(nv,N-1) vs(ns,N) vic(nic) vtc(ntc)

minimize(lx'*(Sx*x(:,N)+cx) + lu'*(Su*u(:,N)+cu) + lp'*(Sp*p+cp) + SCP.lambdaic*norm(vic,2) + SCP.lambdatc*norm(vtc,2) + (norms(Bx*(Sx*x+cx*ones(1,N)),2,1) + norms(Bu*(Su*u+cu*ones(1,N)),2,1) + [SCP.lambda*norms(E(1:nx,1:nv)*v,2,1) 0] + SCP.lambdas*norms(vs,2,1))*[1/2; ones(N-2,1); 1/2]/N)

subject to

% Boundary Conditions:
H0*(Sx*x(1:nx,1)+cx) + K0*(Sp*p+cp) + l0 == 0;
Hf*(Sx*x(1:nx,N)+cx) + Kf*(Sp*p+cp) + lf == 0;

% Convex Constraints on Parameter:
WLp(1:nlp,1:np)*(Sp*p+cp) <= yLp(1:nlp);
for i = 1:length(dqp)
    norm(yQp(sum(dqp(1:i-1))+2:sum(dqp(1:i))) - WQp(sum(dqp(1:i-1))+2:sum(dqp(1:i)),1:np)*(Sp*p+cp),2) <= yQp(sum(dqp(1:i-1))+1) - WQp(sum(dqp(1:i-1))+1,1:np)*(Sp*p+cp);
end

% Trust Region on Parameter:
norm(p-Sp\(pbar-cp),2) <= SCP.etap;

for k = 1:N

    % Dynamics:
    if k < N
        (Sx*x(1:nx,k+1)+cx) == A(1:nx,1:nx,k)*(Sx*x(1:nx,k)+cx) + Bm(1:nx,1:nu,k)*(Su*u(1:nu,k)+cu) + Bp(1:nx,1:nu,k)*(Su*u(1:nu,k+1)+cu) + F(1:nx,1:np,k)*(Sp*p+cp) + r(1:nx,k) + E(1:nx,1:nv)*v(1:nv,k);
    end

    % Convex Constraints
    % Linear Cones:
    WLx(1:nlx,1:nx,k)*(Sx*x(1:nx,k)+cx) <= yLx(1:nlx,k);
    WLu(1:nlu,1:nu,k)*(Su*u(1:nu,k)+cu) <= yLu(1:nlu,k);
    
    % Quadratic Cones:
    for i = 1:length(dqx)
        norm(yQx(sum(dqx(1:i-1))+2:sum(dqx(1:i)),k) - WQx(sum(dqx(1:i-1))+2:sum(dqx(1:i)),1:nx,k)*(Sx*x(1:nx,k)+cx),2) <= yQx(sum(dqx(1:i-1))+1,k) - WQx(sum(dqx(1:i-1))+1,1:nx,k)*(Sx*x(1:nx,k)+cx);
    end
    for i = 1:length(dqu)
        norm(yQu(sum(dqu(1:i-1))+2:sum(dqu(1:i)),k) - WQu(sum(dqu(1:i-1))+2:sum(dqu(1:i)),1:nu,k)*(Su*u(1:nu,k)+cu),2) <= yQu(sum(dqu(1:i-1))+1,k) - WQu(sum(dqu(1:i-1))+1,1:nu,k)*(Su*u(1:nu,k)+cu);
    end

    % Non Convex Constraints:
    C(1:ns,1:nx,k)*(Sx*x(1:nx,k)+cx) + D(1:ns,1:nu,k)*(Su*u(1:nu,k)+cu) + G(1:ns,1:np,k)*(Sp*p+cp) + rp(1:ns,k) + vs(1:ns,k) <= 0;

    % Trust Region on State:
    norm(x(1:nx,k)-Sx\(xbar(1:nx,k)-cx)) + norm(u(1:nu,k)-Su\(ubar(1:nu,k)-cu)) <= SCP.eta;

end

cvx_end

x = Sx*x+cx;
u = Su*u+cu;
p = Sp*p+cp;
L = cvx_optval;

vars.v = v;
vars.vs = vs;
vars.vic = vic;
vars.vtc = vtc;

end

function [Lhat,z] = ConvexSolverMOSEK(c,A,b,G,h,L,Q)

[rcode, res] = mosekopt('symbcon echo(0)');
prob = [];

% Cost:
prob.c = c';

% Linear Inequality:
prob.a = sparse([A;G(1:L,:)]);
prob.buc = [b;h(1:L)];
prob.blc = [b;-inf*ones(L,1)];

% Quadratic Cones:
prob.f = sparse(-G(L+1:end,:));
prob.g = h(L+1:end);
prob.accs = reshape([res.symbcon.MSK_DOMAIN_QUADRATIC_CONE*ones(1,length(Q)); Q'],1,2*length(Q));

% Solve:
[r,res] = mosekopt('minimize',prob);
Lhat = res.sol.itr.pobjval;
z = res.sol.itr.xx;

end

function J = NonLinearCost(SCP,xbar,ubar,pbar,delta,N,lx,lu,lp,Bx,Bu,s,gic,gtc)

% Terminal Cost:
phi = lx'*xbar(:,end) + lu'*ubar(:,end) + lp'*pbar + SCP.lambdaic*norm(gic,2) + SCP.lambdatc*norm(gtc,2);

% Running Cost:
Gamma = (vecnorm(Bx*xbar) + vecnorm(Bu*ubar) + vecnorm(Bh*hbar) + SCP.lambda*vecnorm(delta) + SCP.lambdas*vecnorm(max(s,0))) * [1; ones(N-2,1); 1]*1/N;

% Total Cost:
J = phi + Gamma;

end

function [c,A,b,G,h,L,Q] = Parse(SCP,xbar,ubar,pbar,hbar,lx,lu,lp,Bx,Bu,Bh,A,Bm,Bp,F,H,r,E,WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp,C,D,G,rp,H0,K0,l0,Hf,Kf,lf,Sx,Su,Sp,Sh,cx,cu,cp,ch,N,nx,nu,np,nh,nv,ns,nic,ntc,nlx,nlu,nlp,nlh,nlht,nqx,nqu,nqp,nbx,nbu)

fprintf('Parsing Subproblem... ');

%%%%%% EQUALITY CONSTRAINTS %%%%%%
AM = zeros((N-1)*nx,N*nx); % State Transition Matrix Block
BM = zeros((N-1)*nx,N*nu); % Input Jacobian Block
FM = zeros((N-1)*nx,np); % Parameter Jacobian Block
HM = zeros((N-1)*nx,(N-1)*nh); % Time Varying Parameter Jacobian Block
rv = zeros((N-1)*nx,1); % Residual vector
EM = zeros((N-1)*nx,(N-1)*nv); % Slack Matrix
for i = 1:N-1
    AM(nx*(i-1)+1:nx*i,nx*(i-1)+1:nx*i) = A(1:nx,1:nx,i)*Sx;
    AM(nx*(i-1)+1:nx*i,nx*i+1:nx*(i+1)) = -eye(nx,nx)*Sx;

    BM(nx*(i-1)+1:nx*i,nu*(i-1)+1:nu*i) = Bm(1:nx,1:nu,i)*Su;
    BM(nx*(i-1)+1:nx*i,nu*i+1:nu*(i+1)) = Bp(1:nx,1:nu,i)*Su;

    FM(nx*(i-1)+1:nx*i,1:np) = F(1:nx,1:np,i)*Sp;
    
    HM(nx*(i-1)+1:nx*i,nh*(i-1)+1:nh*i) = H(1:nx,1:nh,i)*Sh;

    rv(nx*(i-1)+1:nx*i) = cx - r(1:nx,i) - A(1:nx,1:nx,i)*cx - Bm(1:nx,1:nu,i)*cu - Bp(1:nx,1:nu,i)*cu - F(1:nx,1:np,i)*cp - H(1:nx,1:nh,i)*cp;

    EM(nx*(i-1)+1:nx*i,nv*(i-1)+1:nv*i) = E(1:nx,1:nv,i);
end
H0M = [H0*Sx zeros(nic,nx*(N-1))]; % Initial Condition Block

HfM = [zeros(ntc,nx*(N-1)) Hf*Sx]; % Final Condition Block

A = [AM BM FM HM EM zeros((N-1)*nx,N*ns) zeros((N-1)*nx,nic) zeros((N-1)*nx,ntc) zeros((N-1)*nx,N) zeros((N-1)*nx,N) zeros((N-1)*nx,N) zeros((N-1)*nx,N) zeros((N-1)*nx,N-1) zeros((N-1)*nx,N) zeros((N-1)*nx,1) zeros((N-1)*nx,1); % Dynamics Equality Constraint
     H0M zeros(nic,N*nu) K0*Sp zeros(nic,(N-1)*nh) zeros(nic,(N-1)*nv) zeros(nic,N*ns) eye(nic,nic) zeros(nic,ntc) zeros(nic,N) zeros(nic,N) zeros(nic,N) zeros(nic,N) zeros(nic,N-1) zeros(nic,N) zeros(nic,1) zeros(nic,1); % Initial Condition Equality Constraint
     HfM zeros(ntc,N*nu) Kf*Sp zeros(nic,(N-1)*nh) zeros(ntc,(N-1)*nv) zeros(ntc,N*ns) zeros(ntc,nic) eye(ntc,ntc) zeros(ntc,N) zeros(ntc,N) zeros(ntc,N) zeros(ntc,N) zeros(ntc,N-1) zeros(ntc,N) zeros(ntc,1) zeros(ntc,1)]; % Final Condition Equality Constraint

b = [rv;-l0-H0*cx-K0*cp;-lf-Hf*cx-Kf*cp];

%%%%%% INEQUALITY CONSTRAINTS %%%%%%

WLxM = zeros(N*nlx,N*nx); % Convex Linear Cone State Block
yLxv = zeros(N*nlx,1); % Convex Linear Cone State Vector
WLuM = zeros(N*nlu,N*nu); % Convex Linear Cone Input Block
yLuv = zeros(N*nlu,1); % Convex Linear Cone Input Vector
WLhM = zeros(N*nlh,N*nh); % Convex Linear Cone Time Varying Parameter Block
yLhv = zeros(N*nlh,1); % Convex Linear Cone Time Varying Parameter Vector
WLhtM = kron(ones(1,(N-1)),WLht*Sh); % Convex Linear Cone Time Varying Parameter Block
yLhtv = yLht-(N-1)*ch; % Convex Linear Cone Time Varying Parameter Vector
CM = zeros(N*ns,N*nx); % Non-Convex State Constraint Block
DM = zeros(N*ns,N*nu); % Non-Convex Input Constraint Block
GM = zeros(N*ns,np); % Non-convex Parameter Constraint Block
IM = zeros(N*ns,nh); % Non-convex Time Varying Parameter Constraint Block
rpv = zeros(N*ns,1); % Non-Convex Residual

WQxM = zeros(N*nqx,N*nx); % Convex Quadratic Cone State Block
yQxv = zeros(N*nqx,1); % Convex Quadratic Cone State Vector
WQuM = zeros(N*nqu,N*nu); % Convex Quadratic Cone Input Block
yQuv = zeros(N*nqu,1); % Convex Quadratic Cone Input Vector

for i = 1:N
    WLxM(nlx*(i-1)+1:nlx*i,nx*(i-1)+1:nx*i) = WLx(1:nlx,1:nx,i)*Sx;
    yLxv(nlx*(i-1)+1:nlx*i) = yLx(1:nlx,i)-WLx(1:nlx,1:nx,i)*cx;
    WLuM(nlu*(i-1)+1:nlu*i,nu*(i-1)+1:nu*i) = WLu(1:nlu,1:nu,i)*Su;
    yLuv(nlu*(i-1)+1:nlu*i) = yLu(1:nlu,i)-WLu(1:nlu,1:nu,i)*cu;
    WLhM(nlh*(i-1)+1:nlh*i,nh*(i-1)+1:nh*i) = WLh(1:nlh,1:nh)*Sh;
    yLhv(nlh*(i-1)+1:nlh*i) = yLh(1:nlh,i)-WLh(1:nlh,1:nh)*ch;

    CM(ns*(i-1)+1:ns*i,nx*(i-1)+1:nx*i) = C(1:ns,1:nx,i)*Sx;
    DM(ns*(i-1)+1:ns*i,nu*(i-1)+1:nu*i) = D(1:ns,1:nu,i)*Su;
    GM(ns*(i-1)+1:ns*i,:) = G(1:ns,1:np,i)*Sp;
    if i < N IM(ns*(i-1)+1:ns*i,nh*(i-1)+1:nh*i) = I(1:ns,1:nh,i)*Sh; end
    rpv(ns*(i-1)+1:ns*i) = -rp(1:ns,i) - C(1:ns,1:nx,i)*cx - D(1:ns,1:nu,i)*cu - G(1:ns,1:np,i)*cp - I(1:ns,1:nh,i)*ch;

    WQxM(nqx*(i-1)+1:nqx*i,nx*(i-1)+1:nx*i) = WQx(1:nqx,1:nx,i)*Sx;
    yQxv(nqx*(i-1)+1:nqx*i) = yQx(1:nqx,i) - WQx(1:nqx,1:nx,i)*cx;
    WQuM(nqu*(i-1)+1:nqu*i,nu*(i-1)+1:nu*i) = WQu(1:nqu,1:nu,i)*Su;
    yQuv(nqu*(i-1)+1:nqu*i) = yQu(1:nqu,i) - WQu(1:nqu,1:nu,i)*cu;
end
yLp = yLp - WLp*cp;
WLp = WLp*Sp;

% State Trust Region Quadratic Cone Constraint
Mxtr = kron(eye(N),[zeros(1,nx); -eye(nx)]);
Nxtr = kron(eye(N),[-1;zeros(nx,1)]);
yxtr = zeros(N*(nx+1),1);
for i = 1:N
    yxtr((i-1)*(nx+1)+1) = 0;
    warning('off');
    yxtr((i-1)*(nx+1)+2:i*(nx+1)) = Sx\(-xbar(:,i)+cx);
    warning('on');
end
if any(isnan(yxtr))
    warning('Scaling Matrices are Singular');
    yxtr(isnan(yxtr)) = 0;
end

% Input Trust Region Quadratic Cone Constraint
Mutr = kron(eye(N),[zeros(1,nu); -eye(nu)]);
Nutr = kron(eye(N),[-1;zeros(nu,1)]);
yutr = zeros(N*(nu+1),1);
for i = 1:N
    yutr((i-1)*(nu+1)+1) = 0;
    warning('off');
    yutr((i-1)*(nu+1)+2:i*(nu+1)) = Su\(-ubar(:,i)+cu);
    warning('on');
end
if any(isnan(yutr))
    warning('Scaling Matrices are Singular');
    yutr(isnan(yutr)) = 0;
end

% Parameter Trust Region Quadratic Cone Constraint
Mptr = [zeros(1,np);-eye(np)];
warning('off');
yptr = [SCP.etap;Sp\(-pbar+cp)];
warning('on');
if any(isnan(yptr))
    warning('Scaling Matrices are Singular');
    yptr(isnan(yptr)) = 0;
end

% State Cost Quadratic Cone Constraint
Mxc = kron(eye(N),[zeros(1,nx);-Bx*Sx]);
Nxc = kron(eye(N),[-1;zeros(nbx,1)]);
yxc = kron(ones(N,1),[0;Bx*cx]);

% Input Cost Quadratic Cone Constraint
Muc = kron(eye(N),[zeros(1,nu);-Bu*Su]);
Nuc = kron(eye(N),[-1;zeros(nbu,1)]);
yuc = kron(ones(N,1),[0;Bu*cu]);

% Dynamic Slack Quadratic Cone Constraint
Mvc = zeros((N-1)*nx,(N-1)*nv);
for i = 1:N-1
    Mvc((i-1)*(nx+1)+1:i*(nx+1),(i-1)*nv+1:i*nv) = [zeros(1,nv);-E(1:nx,1:nv,i)];
end
Nvc = kron(eye(N-1),[-1;zeros(nx,1)]);
yvc = zeros((N-1)*(nx+1),1);

% Non Convex Constraint Slack Quadratic Cone Constraint
Mvsc = kron(eye(N),[zeros(1,ns); -eye(ns)]);
Nvsc = kron(eye(N),[-1;zeros(ns,1)]);
yvsc = zeros(N*(ns+1),1);

% Initial Condition Slack Quadratic Cone Constraint
Mvicc = [zeros(1,nic);-eye(nic)];
Nvicc = [-1; zeros(nic,1)];
yvicc = zeros(nic+1,1);

% Final Condition Slack Quadratic Cone Constraint
Mvtcc = [zeros(1,ntc);-eye(ntc)];
Nvtcc = [-1; zeros(ntc,1)];
yvtcc = zeros(ntc+1,1);

G = [WLxM zeros(N*nlx,N*nu) zeros(N*nlx,np) zeros(N*nlx,(N-1)*nv) zeros(N*nlx,N*ns) zeros(N*nlx,nic) zeros(N*nlx,ntc) zeros(N*nlx,N) zeros(N*nlx,N) zeros(N*nlx,N) zeros(N*nlx,N) zeros(N*nlx,N-1) zeros(N*nlx,N) zeros(N*nlx,1) zeros(N*nlx,1); % Convex State Linear Cone Constraint
     zeros(N*nlu,N*nx) WLuM zeros(N*nlu,np) zeros(N*nlu,(N-1)*nv) zeros(N*nlu,N*ns) zeros(N*nlu,nic) zeros(N*nlu,ntc) zeros(N*nlu,N) zeros(N*nlu,N) zeros(N*nlu,N) zeros(N*nlu,N) zeros(N*nlu,N-1) zeros(N*nlu,N) zeros(N*nlu,1) zeros(N*nlu,1); % Convex Input Linear Cone Constraint
     zeros(nlp,N*nx) zeros(nlp,N*nu) WLp zeros(nlp,(N-1)*nv) zeros(nlp,N*ns) zeros(nlp,nic) zeros(nlp,ntc) zeros(nlp,N) zeros(nlp,N) zeros(nlp,N) zeros(nlp,N) zeros(nlp,N-1) zeros(nlp,N) zeros(nlp,1) zeros(nlp,1); % Convex Parameter Linear Cone Constraint
     CM DM GM zeros(N*ns,(N-1)*nv) -eye(N*ns,N*ns) zeros(N*ns,nic) zeros(N*ns,ntc) zeros(N*ns,N) zeros(N*ns,N) zeros(N*ns,N) zeros(N*ns,N) zeros(N*ns,N-1) zeros(N*ns,N) zeros(N*ns,1) zeros(N*ns,1); % Non-Convex Linear Cone Constraint
     zeros(N,N*nx) zeros(N,N*nu) zeros(N,np) zeros(N,(N-1)*nv) zeros(N,N*ns) zeros(N,nic) zeros(N,ntc) eye(N,N) eye(N,N) zeros(N,N) zeros(N,N) zeros(N,N-1) zeros(N,N) zeros(N,1) zeros(N,1); % Linear Cone State-Input Trust Region Constraint
     WQxM zeros(N*nqx,N*nu) zeros(N*nqx,np) zeros(N*nqx,(N-1)*nv) zeros(N*nqx,N*ns) zeros(N*nqx,nic) zeros(N*nqx,ntc) zeros(N*nqx,N) zeros(N*nqx,N) zeros(N*nqx,N) zeros(N*nqx,N) zeros(N*nqx,N-1) zeros(N*nqx,N) zeros(N*nqx,1) zeros(N*nqx,1); % Convex State Quadratic Cone Constraint
     zeros(N*nqu,N*nx) WQuM zeros(N*nqu,np) zeros(N*nqu,(N-1)*nv) zeros(N*nqu,N*ns) zeros(N*nqu,nic) zeros(N*nqu,ntc) zeros(N*nqu,N) zeros(N*nqu,N) zeros(N*nqu,N) zeros(N*nqu,N) zeros(N*nqu,N-1) zeros(N*nqu,N) zeros(N*nqu,1) zeros(N*nqu,1); % Convex Input Quadratic Cone Constraint
     zeros(nqp,N*nx) zeros(nqp,N*nu) WQp zeros(nqp,(N-1)*nv) zeros(nqp,N*ns) zeros(nqp,nic) zeros(nqp,ntc) zeros(nqp,N) zeros(nqp,N) zeros(nqp,N) zeros(nqp,N) zeros(nqp,N-1) zeros(nqp,N) zeros(nqp,1) zeros(nqp,1); % Convex Parameter Quadratic Cone Constraint
     Mxtr zeros(N*(nx+1),N*nu) zeros(N*(nx+1),np) zeros(N*(nx+1),(N-1)*nv) zeros(N*(nx+1),N*ns) zeros(N*(nx+1),nic) zeros(N*(nx+1),ntc) Nxtr zeros(N*(nx+1),N) zeros(N*(nx+1),N) zeros(N*(nx+1),N) zeros(N*(nx+1),N-1) zeros(N*(nx+1),N) zeros(N*(nx+1),1) zeros(N*(nx+1),1); % State Trust Region Quadratic Cone Constraint
     zeros(N*(nu+1),N*nx) Mutr zeros(N*(nu+1),np) zeros(N*(nu+1),(N-1)*nv) zeros(N*(nu+1),N*ns) zeros(N*(nu+1),nic) zeros(N*(nu+1),ntc) zeros(N*(nu+1),N) Nutr zeros(N*(nu+1),N) zeros(N*(nu+1),N) zeros(N*(nu+1),N-1) zeros(N*(nu+1),N) zeros(N*(nu+1),1) zeros(N*(nu+1),1); % Input Trust Region Quadratic Cone Constraint
     zeros(np+1,N*nx) zeros(np+1,N*nu) Mptr zeros(np+1,(N-1)*nv) zeros(np+1,N*ns) zeros(np+1,nic) zeros(np+1,ntc) zeros(np+1,N) zeros(np+1,N) zeros(np+1,N) zeros(np+1,N) zeros(np+1,N-1) zeros(np+1,N) zeros(np+1,1) zeros(np+1,1); % Parameter Trust Region Quadratic Cone Constraint
     Mxc zeros(N*(nbx+1),N*nu) zeros(N*(nbx+1),np) zeros(N*(nbx+1),(N-1)*nv) zeros(N*(nbx+1),N*ns) zeros(N*(nbx+1),nic) zeros(N*(nbx+1),ntc) zeros(N*(nbx+1),N) zeros(N*(nbx+1),N) Nxc zeros(N*(nbx+1),N) zeros(N*(nbx+1),N-1) zeros(N*(nbx+1),N) zeros(N*(nbx+1),1) zeros(N*(nbx+1),1); % State Cost Quadratic Cone Constraint
     zeros(N*(nbu+1),N*nx) Muc zeros(N*(nbu+1),np) zeros(N*(nbu+1),(N-1)*nv) zeros(N*(nbu+1),N*ns) zeros(N*(nbu+1),nic) zeros(N*(nbu+1),ntc) zeros(N*(nbu+1),N) zeros(N*(nbu+1),N) zeros(N*(nbu+1),N) Nuc zeros(N*(nbu+1),N-1) zeros(N*(nbu+1),N) zeros(N*(nbu+1),1) zeros(N*(nbu+1),1); % Input Cost Quadratic Cone Constraint
     zeros((N-1)*(nx+1),N*nx) zeros((N-1)*(nx+1),N*nu) zeros((N-1)*(nx+1),np) Mvc zeros((N-1)*(nx+1),N*ns) zeros((N-1)*(nx+1),nic) zeros((N-1)*(nx+1),ntc) zeros((N-1)*(nx+1),N) zeros((N-1)*(nx+1),N) zeros((N-1)*(nx+1),N) zeros((N-1)*(nx+1),N) Nvc zeros((N-1)*(nx+1),N) zeros((N-1)*(nx+1),1) zeros((N-1)*(nx+1),1); % Dynamic Slack Cost Quadratic Cone Constraint
     zeros(N*(ns+1),N*nx) zeros(N*(ns+1),N*nu) zeros(N*(ns+1),np) zeros(N*(ns+1),(N-1)*nv) Mvsc zeros(N*(ns+1),nic) zeros(N*(ns+1),ntc) zeros(N*(ns+1),N) zeros(N*(ns+1),N) zeros(N*(ns+1),N) zeros(N*(ns+1),N) zeros(N*(ns+1),N-1) Nvsc zeros(N*(ns+1),1) zeros(N*(ns+1),1); % Non-Convex Slack Cost Quadratic Cone Constraint
     zeros(nic+1,N*nx) zeros(nic+1,N*nu) zeros(nic+1,np) zeros(nic+1,(N-1)*nv) zeros(nic+1,N*ns) Mvicc zeros(nic+1,ntc) zeros(nic+1,N) zeros(nic+1,N) zeros(nic+1,N) zeros(nic+1,N) zeros(nic+1,N-1) zeros(nic+1,N) Nvicc zeros(nic+1,1); % Initial Condition Slack Cost Quadratic Cone Constraint
     zeros(ntc+1,N*nx) zeros(ntc+1,N*nu) zeros(ntc+1,np) zeros(ntc+1,(N-1)*nv) zeros(ntc+1,N*ns) zeros(ntc+1,nic) Mvtcc zeros(ntc+1,N) zeros(ntc+1,N) zeros(ntc+1,N) zeros(ntc+1,N) zeros(ntc+1,N-1) zeros(ntc+1,N) zeros(ntc+1,1) Nvtcc; % Initial Condition Slack Cost Quadratic Cone Constraint
    ];

h = [yLxv;yLuv;yLp;rpv;ones(N,1)*SCP.eta;yQxv;yQuv;yQp;yxtr;yutr;yptr;yxc;yuc;yvc;yvsc;yvicc;yvtcc];

% Cone Dimensions
L = N*(nlx+nlu+ns+1) + nlp; % Linear cone Dimensions
Q = [kron(ones(N,1),dqx); kron(ones(N,1),dqu); dqp; (nx+1)*ones(N,1); (nu+1)*ones(N,1); np+1; (nbx+1)*ones(N,1); (nbu+1)*ones(N,1); (nx+1)*ones(N-1,1); (ns+1)*ones(N,1); nic+1; ntc+1]; % Quadratic Cone Dimensions


%%%%%% AFFINE COST FUNCTION %%%%%%
c = [[zeros((N-1)*nx,1); (lx'*Sx)']; % State
     [zeros((N-1)*nu,1); (lu'*Su)']; % Input
     (lp'*Sp)'; % Parameter
     zeros((N-1)*nv,1); % Dynamic Slack
     zeros(N*ns,1); % Constraint Slack
     zeros(nic,1); % Initial Condition Slack
     zeros(ntc,1); % Final Condition Slack
     zeros(N,1); % State Trust Region
     zeros(N,1); % Input Trust Region
     1/N*[1; ones(N-2,1); 1]; % State Integral
     1/N*[1; ones(N-2,1); 1]; % Input Integral
     SCP.lambda*1/N*[1; ones(N-3,1); 1]; % Dynamic Slack Penalty
     SCP.lambdas*1/N*[1; ones(N-2,1); 1]; % Constraint Slack Penalty
     SCP.lambdaic; % Initial Condition Slack
     SCP.lambdatc]; % Final Condition Slack

fprintf('Done\n');

end

function [L,x,u,p,vars] = Descale(Lhat,z,Sx,cx,Su,cu,Sp,cp,lx,lu,lp,N,nx,nu,np,nv,ns,nic,ntc)

% Decompose Solution
x = Sx*reshape(z(1:nx*N),nx,N)+cx; % State
u = Su*reshape(z(nx*N+1:nx*N+nu*N),nu,N)+cu; % Input
p = Sp*z(nx*N+nu*N+1:nx*N+nu*N+np)+cp; % Parameter
L = Lhat + lx'*cx+lu'*cu+lp'*cp; % Linear Cost

% Additional Solution Variables:
vars.v = reshape(z(nx*N+nu*N+np+1:nx*N+nu*N+np+nv*(N-1)),nv,N-1);
vars.vs = reshape(z(nx*N+nu*N+np+nv*(N-1)+1:nx*N+nu*N+np+nv*(N-1)+ns*N),ns,N);
vars.vic = z(nx*N+nu*N+np+nv*(N-1)+ns*N+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic);
vars.vtc = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc);
vars.deltax = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+N);
vars.deltau = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+N+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+2*N);
vars.bxnorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+2*N+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+3*N);
vars.bunorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+3*N+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+4*N);
vars.evnorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+4*N+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+4*N+(N-1));
vars.vsnorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+4*N+(N-1)+1:nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+5*N+(N-1));
vars.vicnorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+5*N+(N-1)+1);
vars.vtcnorm = z(nx*N+nu*N+np+nv*(N-1)+ns*N+nic+ntc+5*N+(N-1)+2);

end

function [A,Bm,Bp,F,r,E,delta] = DiscretizeDynamics(SCP,System,xbar,ubar,pbar,t,nx,nu,np,nh,nv,N,Dynamics)

fprintf('Discretizing Subproblem... ')

% Preallocate:
delta = zeros(nx,N);
A = zeros(nx,nx,N);
Bm = zeros(nx,nu,N);
Bp = zeros(nx,nu,N);
F = zeros(nx,np,N);
H = zeros(nx,nh,N);
r = zeros(nx,N);
E = zeros(nx,nv,N);

% For Loop to Compute Discrete Time Dynamics at Each Time Step
for k = 1:N-1
    
    % Initial Conditions
    phi0 = eye(nx); % Initial Condition for A
    PBm0 = zeros(nx,nu); % Initial Condition for Bm
    PBp0 = zeros(nx,nu); % Initial Condition for Bp
    PF0 = zeros(nx,np);  % Initial Condition for F
    PH0 = zeros(nx,nh); % Initial COndition for H
    Pr0 = zeros(nx,1); % Initial Condition for r

    % Vectorize Initial Conditions
    P0 = Flatten(xbar(:,k),phi0,PBm0,PBp0,PF0,PH0,Pr0,nx,nu,np,nh);
    
    % Integrate from current time step to next:
    tspan = linspace(t(k),t(k+1),SCP.Nsub);
    P = RK4(@(tc,x) Derivatives(System,x,tc,ubar(:,k+1),ubar(:,k),pbar,hbar,t(k+1),t(k),nx,nu,np,nh,Dynamics),tspan,P0);
    
    % Unvectorize Integrated Result
    [xf,Pphi,PBm,PBp,Ps,Pr] = Unflatten(P,nx,nu,np);
    
    % Discretization Defect
    delta(1:nx,k) = xbar(:,k+1)-xf;
    
    % Compute Discrete Time Dynamics Matrices
    A(1:nx,1:nx,k) = Pphi;
    Bm(1:nx,1:nu,k) = Pphi*PBm;
    Bp(1:nx,1:nu,k) = Pphi*PBp;
    F(1:nx,1:np,k) = Pphi*Ps;
    r(1:nx,k) = Pphi*Pr;
    [~,~,~,~,E(1:nx,1:nv,k)] = Dynamics(System,xbar(:,k),ubar(:,k),pbar,t(k));

end

fprintf('Done \n');

end

function [x,Pphi,PBm,PBp,PF,PH,Pr] = Unflatten(P,nx,nu,np,nh)

x = P(1:nx);
Pphi = reshape(P(nx+1:nx+nx^2),[nx,nx]);
PBm = reshape(P(nx+nx^2+1:nx+nx^2+nx*nu),[nx,nu]);
PBp = reshape(P(nx+nx^2+nx*nu+1:nx+nx^2+2*nx*nu),[nx,nu]);
PF = reshape(P(nx+nx^2+2*nx*nu+1:nx+nx^2+2*nx*nu+nx*np),[nx,np]);
PH = reshape(P(nx+nx^2+2*nx*nu+nx*np+1:nx+nx^2+2*nx*nu+nx*np+nx*nh),[nx,nh]);
Pr = P(nx+nx^2+2*nx*nu+nx*np+nx*nh:2*nx+nx^2+2*nx*nu+nx*np+nx*nh);

end

function P = Flatten(x,Phi,PBm,PBp,PF,Pr,nx,nu,np,nh)

P = [x;
     reshape(Phi,[nx^2,1]);
     reshape(PBm,[nx*nu,1]);
     reshape(PBp,[nx*nu,1]);
     reshape(PF,[nx*np,1]);
     reshape(PH,[nx*nh,1]);
     Pr];

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

function [Pdot] = Derivatives(System,P,t,up,um,p,tp,tm,nx,nu,np,Dynamics)

% Unvectorize Integration Vector
[x,Pphi,~,~,~,~,~] = Unflatten(P,nx,nu,np,nh);

% Compute Time Ratios in Integration for First Order Hold
lambdakm = (tp-t)/(tp-tm);
lambdakp = 1-lambdakm;

% Determine Input from First Order Hold
u = lambdakm*um + lambdakp*up;

% Compute Dynamics
[f,A,B,F,H,~] = Dynamics(System,x,u,p,h,t);
r = f - A*x - B*u - F*p - H*h; % Residual

% Invert State Transition Matrix
psi = Pphi^-1;

% Compute Derivatives 
Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PFdot = psi*F;
PHdot = psi*H;
PRdot = psi*r;

% Vectorize Result
Pdot = Flatten(Pxdot,Pphidot,PBmdot,PBpdot,PFdot,PHdot,PRdot,nx,nu,np);

end

function [WLx,yLx,WLu,yLu,WLp,yLp,WLh,yLh,WLht,yLht,WQx,yQx,WQu,yQu,WQp,yQp,dqx,dqu,dqp] = DiscretizeConvexConstraint(System,xbar,ubar,pbar,t,nx,nu,nlx,nlu,nlh,nlht,nqx,nqu,N,ConvexConstraints)

% Preallocate
WLx = zeros(nlx,nx,N);
yLx = zeros(nlx,N);
WLu = zeros(nlu,nu,N);
yLu = zeros(nlu,N);
WLh = zeros(nlh,nh,N);
yLh = zeros(nlh,N);
WLht = zeros(nlht,nh,N);
yLht = zeros(nlht,N);
WQx = zeros(nqx,nx,N);
yQx = zeros(nqx,N);
WQu = zeros(nqu,nu,N);
yQu = zeros(nqu,N);

% Discretize
for i = 1:N

    [WLx(1:nlx,1:nx,i),yLx(1:nlx,i),WLu(1:nlu,1:nu,i),yLu(1:nlu,i),WLp,yLp,WLh(1:nlh,nh,i),yLh(1:nlh,i),WLht(1:nlht,1:nh,i),yLht(1:nlht,i),WQx(1:nqx,1:nx,i),yQx(1:nqx,i),WQu(1:nqu,1:nu,i),yQu(1:nqu,i),WQp,yQp,dqx,dqu,dqp] = ConvexConstraints(System,xbar(:,i),ubar(:,i),pbar,hbar(:,i),t(i));

end

end

function [s,C,D,G,I,rp] = DiscretizeNonConvexConstraints(System,xbar,ubar,pbar,hbar,t,nx,nu,np,nh,ns,N,NonConvexConstraints)

% Preallocate
s = zeros(ns,N);
C = zeros(ns,nx,N);
D = zeros(ns,nu,N);
G = zeros(ns,np,N);
I = zeros(ns,nh,N);
rp = zeros(ns,N);

% Discretize
for i = 1:N

    [s(1:ns,i),C(1:ns,1:nx,i),D(1:ns,1:nu,i),G(1:ns,1:np,i),I(1:ns,1:nh,i)] = NonConvexConstraints(System,xbar(:,i),ubar(:,i),pbar,hbar(:,i),t(i));
    rp(1:ns,i) = s(1:ns,i) - C(1:ns,1:nx,i)*xbar(:,i) - D(1:ns,1:nu,i)*ubar(:,i) - G(1:ns,1:np,i)*pbar - I(1:ns,1:nh,i)*hbar(:,i);

end

end

function [Sx,Su,Sp,Sh,cx,cu,cp,ch] = ScalingMatrices(System,Scaling)

[xmin,xmax,umin,umax,pmin,pmax] = Scaling(System);

Sx = diag(xmax-xmin);
cx = xmin;

Su = diag(umax-umin);
cu = umin;

Sp = diag(pmax-pmin);
cp = pmin;

Sh = diag(hmax-hmin);
ch = hmin;

end