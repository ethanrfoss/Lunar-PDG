function PlotRocket(System,x,u)

% Dimensions
cent = 100/1000; %Centroid[m]
d = 20/1000; %Body Diameter[m]
ds = 10/1000; %Skirt Diameter[m]
lc = 40/1000; %Cone Length[m]
lb = 100/1000; %Body Length[m]
ls = 20/1000; %Skirt Length[m]
lm = 10/1000; %Engine Length[m]
lf = 40/1000; %Fire Length Max[m]

%Cone
cone.X = [0 d/2*cos(0:pi/16:2*pi)];
cone.Y = [0 d/2*sin(0:pi/16:2*pi)];
cone.Z = [0 -lc*ones(1,length(0:pi/16:2*pi))];
cone.hull = convhull(cone.X,cone.Y,cone.Z);

%Body
body.X = [d/2*cos(0:pi/16:2*pi) d/2*cos(0:pi/16:2*pi)];
body.Y = [d/2*sin(0:pi/16:2*pi) d/2*sin(0:pi/16:2*pi)];
body.Z = [-lc*ones(1,length(0:pi/16:2*pi)) -(lc+lb)*ones(1,length(0:pi/16:2*pi))];
body.hull = convhull(body.X,body.Y,body.Z);

%Skirt
skirt.X = [d/2*cos(0:pi/16:2*pi) ds/2*cos(0:pi/16:2*pi)];
skirt.Y = [d/2*sin(0:pi/16:2*pi) ds/2*sin(0:pi/16:2*pi)];
skirt.Z = [-(lc+lb)*ones(1,length(0:pi/16:2*pi)) -(lc+lb+ls)*ones(1,length(0:pi/16:2*pi))];
skirt.hull = convhull(skirt.X,skirt.Y,skirt.Z);

%Motor
motor.X = [ds/2*cos(0:pi/16:2*pi) 0];
motor.Y = [ds/2*sin(0:pi/16:2*pi) 0];
motor.Z = [-(lc+lb+ls+lm)*ones(1,length(0:pi/16:2*pi)) -(lc+lb+ls)];
motor.hull = convhull(motor.X,motor.Y,motor.Z);

%Fire
fire.X = [ds/2*cos(0:pi/16:2*pi) 0];
fire.Y = [ds/2*sin(0:pi/16:2*pi) 0];
fire.Z = [-(lc+lb+ls+lm)*ones(1,length(0:pi/16:2*pi)) -(lc+lb+ls+lm+lf)];
fire.hull = convhull(fire.X,fire.Y,fire.Z);

Scale = 100;
AScale = 1;

% Rotation:
TB2L = @(q) [1-2*(q(3)^2+q(4)^2) 2*(q(2)*q(3)+q(4)*q(1)) 2*(q(2)*q(4)-q(3)*q(1));
             2*(q(2)*q(3)-q(4)*q(1)) 1-2*(q(2)^2+q(4)^2) 2*(q(3)*q(4)+q(2)*q(1));
             2*(q(2)*q(4)+q(3)*q(1)) 2*(q(3)*q(4)-q(2)*q(1)) 1-2*(q(2)^2+q(3)^2)]';

% Plot:
figure; hold on; axis equal; grid on; axis tight;
plot3(x(1,:),x(2,:),x(3,:),'Color','k','LineWidth',2);
for i = 1:length(x(1,:))
    
    xb = cross([1;0;0],u(1:3,i))/norm(cross([1;0;0],u(1:3,i)));
    Ta = [xb cross(xb,u(1:3,i))/norm(cross(xb,u(1:3,i))) u(1:3,i)/norm(u(1:3,i))];
    TBL = TB2L(x(7:10,i));

    Motor = Ta*([motor.X;motor.Y;motor.Z]-[0;0;-(lc+lb+ls)])+[0;0;-(lc+lb+ls)];
    Fire = Ta*([fire.X;fire.Y;[fire.Z(1:end-1) -(lc+lb+ls+lm+lf*u(4,i)/System.MaximumThrust)]]-[0;0;-(lc+lb+ls)])+[0;0;-(lc+lb+ls)];
    Cone = Scale*TBL*([cone.X;cone.Y;cone.Z]-[0;0;-cent])+[x(1,i);x(2,i);x(3,i)];
    Body = Scale*TBL*([body.X;body.Y;body.Z]-[0;0;-cent])+[x(1,i);x(2,i);x(3,i)];
    Skirt = Scale*TBL*([skirt.X;skirt.Y;skirt.Z]-[0;0;-cent])+[x(1,i);x(2,i);x(3,i)];
    Motor = Scale*TBL*(Motor-[0;0;-cent])+[x(1,i);x(2,i);x(3,i)];
    Fire = Scale*TBL*(Fire-[0;0;-cent])+[x(1,i);x(2,i);x(3,i)];

    trisurf(cone.hull,Cone(1,:),Cone(2,:),Cone(3,:),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
    trisurf(body.hull,Body(1,:),Body(2,:),Body(3,:),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
    trisurf(skirt.hull,Skirt(1,:),Skirt(2,:),Skirt(3,:),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
    trisurf(motor.hull,Motor(1,:),Motor(2,:),Motor(3,:),'FaceColor','k','EdgeColor','k');
    trisurf(fire.hull,Fire(1,:),Fire(2,:),Fire(3,:),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
    plot3(Cone(1,:),Cone(2,:),Cone(3,:),'k');
    plot3(Body(1,:),Body(2,:),Body(3,:),'k');
    plot3(Skirt(1,:),Skirt(2,:),Skirt(3,:),'k');

end

axis tight;

end