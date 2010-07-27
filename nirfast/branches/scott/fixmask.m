function pts = fixmask(mesh)

plotimage(mesh,mesh.region);
[x,y] = ginput;

for i = 1:length(x);
    pts(i) = find(((x(i)-mesh.nodes(:,1)).^2+(y(i)-mesh.nodes(:,2)).^2).^0.5 ==(min(((x(i)-mesh.nodes(:,1)).^2+(y(i)-mesh.nodes(:,2)).^2).^0.5)));
end