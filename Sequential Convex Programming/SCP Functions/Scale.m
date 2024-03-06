function [ScalingMatrices] = Scale(SCP,System,Scaling)

[xmin,xmax,umin,umax] = Scaling(System);

ScalingMatrices.Sx = diag(xmax-xmin) + diag(xmax==xmin);
ScalingMatrices.cx = xmin;

ScalingMatrices.Su = diag(umax-umin) + diag(umax==umin);
ScalingMatrices.cu = umin;

if SCP.AdaptiveMesh
    ScalingMatrices.Sh = System.MaxStep - System.MinStep + (System.MaxStep == System.MinStep);
    ScalingMatrices.ch = System.MinStep;
end
if SCP.FreeTime && ~SCP.AdaptiveMesh
    ScalingMatrices.Sp = 1;%System.tmax-System.tmin + (System.tmax==System.tmin);
    ScalingMatrices.cp = System.tmin;
end

end