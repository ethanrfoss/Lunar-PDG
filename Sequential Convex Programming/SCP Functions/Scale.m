function [ScalingMatrices] = Scale(System,Scaling)

[xmin,xmax,umin,umax,pmin,pmax,hmin,hmax] = Scaling(System);

ScalingMatrices.Sx = diag(xmax-xmin) + diag(xmax==xmin);
ScalingMatrices.cx = xmin;

ScalingMatrices.Su = diag(umax-umin) + diag(umax==umin);
ScalingMatrices.cu = umin;

ScalingMatrices.Sp = diag(pmax-pmin) + diag(pmax==pmin);
ScalingMatrices.cp = pmin;

ScalingMatrices.Sh = diag(hmax-hmin) + diag(hmax == hmin);
ScalingMatrices.ch = hmin;

end