
function [Result] = ConvexSolve(ConvexProblem)

% Solve Cone Program

[r,Result] = mosekopt('minimize echo(0)',ConvexProblem);

if strcmp(Result.sol.itr.prosta,'PRIMAL_AND_DUAL_INFEASIBLE')
    error("Primal and Dual Problem Infeasible");
elseif strcmp(Result.sol.itr.prosta,'PRIMAL_INFEASIBLE')
    error("Primal Problem Infeasible");
elseif strcmp(Result.sol.itr.prosta,'DUAL_INFEASIBLE')
    error("Dual Problem Infeasible");
else
    if ~strcmp(Result.sol.itr.solsta,'OPTIMAL')
        warning('Convex Subproblem is Non-optimal');
    end
end

end