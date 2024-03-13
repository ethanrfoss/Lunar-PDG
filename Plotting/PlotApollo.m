function PlotApollo(System,x,u)

% Dimensions
hb11 = 3;
lb11 = 4.5;
wb11 = 2.75;
hb12 = 1.5;
lb12 = 2;
wb12 = 4.5;
hb2 = 2;
lb2 = 2;
vl1 = [3 4.5;
       0 0;
       0 -3.5];
vl2 = [-3 -4.5;
       0 0;
       0 -3.5];
vl3 = [0 0;
       3 4.5;
       0 -3.5];
vl4 = [0 0;
       -3 -4.5;
       0 -3.5];
tl(:,1) = [2.5;2.5;2];
td(:,1,1) = [1;0;0];
td(:,2,1) = [0;1;0];
td(:,3,1) = [0;0;1];
tl(:,2) = [-2.5;2.5;2];
td(:,1,2) = [-1;0;0];
td(:,2,2) = [0;1;0];
td(:,3,2) = [0;0;1];
tl(:,3) = [2.5;-2.5;2];
td(:,1,3) = [1;0;0];
td(:,2,3) = [0;-1;0];
td(:,3,3) = [0;0;1];
tl(:,4) = [-2.5;-2.5;2];
td(:,1,4) = [-1;0;0];
td(:,2,4) = [0;-1;0];
td(:,3,4) = [0;0;1];
rrcs = .25;
lrcs = .6;
lrcsf = 1;
rc = 1.5;
lc = 2.5;
lf = 3;

% Body 11
Body11.X = [-lb11 lb11 lb11 -lb11 -lb11 -lb11 lb11 lb11 -lb11 -lb11]/2;
Body11.Y = [-wb11 -wb11 wb11 wb11 -wb11 -wb11 -wb11 wb11 wb11 -wb11]/2;
Body11.Z = [0*ones(1,5) hb11*ones(1,5)];
Body11.hull = convhull(Body11.X,Body11.Y,Body11.Z);

% Body 12
Body12.X = [-lb12 lb12 lb12 -lb12 -lb12 -lb12 lb12 lb12 -lb12 -lb12]/2;
Body12.Y = [-wb12 -wb12 wb12 wb12 -wb12 -wb12 -wb12 wb12 wb12 -wb12]/2;
Body12.Z = [0*ones(1,5) hb12*ones(1,5)];
Body12.hull = convhull(Body12.X,Body12.Y,Body12.Z);

% Body 2
Body2.X = [cos(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),sin(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),cos(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),sin(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8))];
Body2.Y = [sin(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),cos(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),sin(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8)),cos(linspace(pi/8,2*pi+pi/8,9))*lb2/(2*sin(pi/8))];
Body2.Z = [0*ones(1,18) -hb2*ones(1,18)];
Body2.hull = convhull(Body2.X,Body2.Y,Body2.Z);

% DPSThruster
DPSThruster = Cone([0;0;0],[0;0;-1],rc,lc);

% RCSThrusters
for i = 1:4
    for j = 1:3
        RCSThruster(i,j) = Cone(tl(:,i),td(:,j,i),rrcs,lrcs);
    end
end

Scale = 2;

% Rotation:
TB2L = @(q) [1-2*(q(3)^2+q(4)^2) 2*(q(2)*q(3)+q(4)*q(1)) 2*(q(2)*q(4)-q(3)*q(1));
             2*(q(2)*q(3)-q(4)*q(1)) 1-2*(q(2)^2+q(4)^2) 2*(q(3)*q(4)+q(2)*q(1));
             2*(q(2)*q(4)+q(3)*q(1)) 2*(q(3)*q(4)-q(2)*q(1)) 1-2*(q(2)^2+q(3)^2)]';

% Plot:
figure; hold on; axis equal; grid on; axis tight;
plot3(x(1,:),x(2,:),x(3,:),'Color','k','LineWidth',2);
for i = 1:length(x(1,:))
    
    TBL = TB2L(x(7:10,i));
    
    lfire = min(max(u(1,i)/System.DPSMaximumThrust,.001),1)*lf;
    DPSFire = Cone([0;0;-lc-lfire],[0;0;1],rc,lfire);
    DPSF = Scale*TBL*[DPSFire.X;DPSFire.Y;DPSFire.Z] + [x(1,i);x(2,i);x(3,i)];
    trisurf(DPSFire.hull,DPSF(1,:),DPSF(2,:),DPSF(3,:),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
    for j = 1:4
        for k = 1:3
            lrcsfire = min(max(u(1+j*k,i)/System.RCSMaximumThrust,.001),1)*lrcsf;
            RCSFire(j,k) = Cone(tl(:,j)+td(:,k,j)*(lrcs+lrcsfire),-td(:,k,j),rrcs,lrcsfire);
            RCSF = Scale*TBL*[RCSFire(j,k).X;RCSFire(j,k).Y;RCSFire(j,k).Z] + [x(1,i);x(2,i);x(3,i)];
            trisurf(RCSFire(j,k).hull,RCSF(1,:),RCSF(2,:),RCSF(3,:),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
            RCS = Scale*TBL*[RCSThruster(j,k).X;RCSThruster(j,k).Y;RCSThruster(j,k).Z] + [x(1,i);x(2,i);x(3,i)];
            trisurf(RCSThruster(j,k).hull,RCS(1,:),RCS(2,:),RCS(3,:),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
        end
    end
    
    DPS = Scale*TBL*([DPSThruster.X;DPSThruster.Y;DPSThruster.Z]) + [x(1,i);x(2,i);x(3,i)];
    trisurf(DPSThruster.hull,DPS(1,:),DPS(2,:),DPS(3,:),'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
    B11 = Scale*TBL*([Body11.X;Body11.Y;Body11.Z]) + [x(1,i);x(2,i);x(3,i)];
    trisurf(Body11.hull,B11(1,:),B11(2,:),B11(3,:),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
    B12 = Scale*TBL*([Body12.X;Body12.Y;Body12.Z]) + [x(1,i);x(2,i);x(3,i)];
    trisurf(Body12.hull,B12(1,:),B12(2,:),B12(3,:),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
    B2 = Scale*TBL*([Body2.X;Body2.Y;Body2.Z]) + [x(1,i);x(2,i);x(3,i)];
    trisurf(Body2.hull,B2(1,:),B2(2,:),B2(3,:),'FaceColor',[1 .843 0],'EdgeColor',[0 0 0]);
    vl1r = Scale*TBL*vl1 + [x(1,i);x(2,i);x(3,i)];
    plot3(vl1r(1,:),vl1r(2,:),vl1r(3,:),'LineWidth',2,'Color',[0 0 0]);
    vl2r = Scale*TBL*vl2 + [x(1,i);x(2,i);x(3,i)];
    plot3(vl2r(1,:),vl2r(2,:),vl2r(3,:),'LineWidth',2,'Color',[0 0 0]);
    vl3r = Scale*TBL*vl3 + [x(1,i);x(2,i);x(3,i)];
    plot3(vl3r(1,:),vl3r(2,:),vl3r(3,:),'LineWidth',2,'Color',[0 0 0]);
    vl4r = Scale*TBL*vl4 + [x(1,i);x(2,i);x(3,i)];
    plot3(vl4r(1,:),vl4r(2,:),vl4r(3,:),'LineWidth',2,'Color',[0 0 0]);
    
    plot3(B11(1,:),B11(2,:),B11(3,:),'k','LineWidth',2);
    plot3(B12(1,:),B12(2,:),B12(3,:),'k','LineWidth',2);
    plot3(B2(1,:),B2(2,:),B2(3,:),'k','LineWidth',2);

end

axis tight;

end

function cone = Cone(loc,dir,rad,len)

T = zeros(3,3);
T(:,3) = dir/norm(dir);
t = cross([sqrt(2)/2;sqrt(2)/2;0],dir)/norm(cross([sqrt(2)/2;sqrt(2)/2;0],dir));
T(:,1) = t;
T(:,2) = cross(T(:,3),T(:,1));

C = [0 rad*cos(linspace(0,2*pi,100));0 rad*sin(linspace(0,2*pi,100));0 len*ones(1,100)];
C = T*C;

cone.X = C(1,:)+loc(1);
cone.Y = C(2,:)+loc(2);
cone.Z = C(3,:)+loc(3);
cone.hull = convhull(cone.X,cone.Y,cone.Z);

end