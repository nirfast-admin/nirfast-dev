function plot_nodes(mesh)


figure;
ind = find(mesh.bndvtx==1);
plot(mesh.source.coord(:,1),mesh.source.coord(:,2),'ro','LineWidth',2,'MarkerSize',8);
hold on;
plot(mesh.meas.coord(:,1),mesh.meas.coord(:,2),'bx','LineWidth',2,'MarkerSize',8);
plot(mesh.nodes(ind,1),mesh.nodes(ind,2),'c.');
axis equal;
legend('Source','Detector');
trisurf(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3))
view(2)
axis equal
