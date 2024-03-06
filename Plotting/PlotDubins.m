
function PlotDubins(System,Boundary,x)

figure; hold on; axis equal; axis tight;
plot(polyshape(Boundary(1,:),Boundary(2,:)),'FaceAlpha',0,'EdgeColor','k');
plot(x(1,:),x(2,:));
Circle = polyshape(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)));
for i = 1:length(System.Obstacles)
    plot(polyshape(System.Obstacles(i).r*Circle.Vertices+System.Obstacles(i).p'),'EdgeColor','k','FaceColor','b');
end
for i = 1:length(x(1,:))
    plot(polyshape(System.CarRadius*Circle.Vertices+x(1:2,i)'),'FaceAlpha',0,'EdgeColor','k');
    plot([x(1,i),x(1,i)+cos(x(3,i))*System.CarRadius],[x(2,i),x(2,i,end)+sin(x(3,i))*System.CarRadius],'k');
end

end