% Places sources and detectors for spec/tomo system.
% Must have SD txt file to load with 4 points at the locations of the
% fiducial markers from the fiber holder.
% M. Mastanduno 6/20/10


% load fiducial points
sd_fn = 'gel_corner.txt';
mesh_fn = 'gel_anom';

mesh = load_mesh(mesh_fn);


fid=fopen(sd_fn,'r');
for i=1:12
    tline=fgetl(fid);
end
a=fgetl(fid);
x= [];y=[];z=[];
while(a~=-1)
    str = str2num(a(6:end));
    x(int32(str(1))) = str(2);
    y(int32(str(1))) = str(3);
    z(int32(str(1))) = str(4);
    a = fgetl(fid);
end
x = x';
y = y';
z = z';
fclose(fid);

if length(x)>4
    cor_mim = [x(5:8) y(5:8) z(5:8)];
    f = [x(1:4) y(1:4) z(1:4)];
else
    f = [x y z];
end

if x(1)==x(2)
    f = [f(:,2) f(:,3)];
elseif y(1)==y(2)
    f = [f(:,1) f(:,3)];
else % z(1)==z(2) 
    f = [f(:,1) f(:,2)];
end

if length(x)>4
    if x(1)==x(2)
        cor_mim = [cor_mim(:,2) cor_mim(:,3)];
    elseif y(1)==y(2)
        cor_mim = [cor_mim(:,1) cor_mim(:,3)];
    else % z(1)==z(2)
        cor_mim = [cor_mim(:,1) cor_mim(:,2)];
    end
end

% account for bitmap transformation
[i j] = max(mesh.nodes(:,2)); % upper left
cor_bmp(1,:) = mesh.nodes(j,1:2);
[i j] = min(mesh.nodes(:,1)); % lower left
cor_bmp(2,:) = mesh.nodes(j,1:2);
[i j] = min(mesh.nodes(:,2)); % lower right
cor_bmp(3,:) = mesh.nodes(j,1:2);
[i j] = max(mesh.nodes(:,1)); % upper right
cor_bmp(4,:) = mesh.nodes(j,1:2);

trans = cor_mim-cor_bmp;
trans = mean(trans);
fid_pts = f-repmat(trans,4,1);

clf
trimesh(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2))
hold on;
plot(cor_mim(:,1),cor_mim(:,2),'ro','linewidth',2)
plot(cor_bmp(:,1),cor_bmp(:,2),'bo','linewidth',2)
plot(fid_pts(:,1),fid_pts(:,2),'go','linewidth',2)
% plot(f(:,1),f(:,2),'ro')
% plot(fid_pts(:,1),fid_pts(:,2),'bo')

dist = 8.9 + (0:12.7:7*12.7);

% place 1-8
shift(1:8,:) = repmat(f(1,:),8,1);
move1 = f(2,:) - f(1,:);
[theta1 rho1] = cart2pol(move1(1),move1(2));
sd(1:8,1) = theta1;
sd(1:8,2) = dist;
% place 9-16
shift(9:16,:) = repmat(f(4,:),8,1);
move2 = f(3,:) - f(4,:);
[theta2,rho2] = cart2pol(move2(1),move2(2));
sd(16:-1:9,1) = theta2;
sd(16:-1:9,2) = dist;
% back to cartesian
[sd(:,1) sd(:,2)] = pol2cart(sd(:,1),sd(:,2));
sd = sd+shift;
% dicom to bmp transform
sd = sd-repmat(trans,16,1);


clf
trimesh(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2))
hold on;
plot(sd(:,1),sd(:,2),'ro','linewidth',2)
plot(fid_pts(:,1),fid_pts(:,2),'go','linewidth',2)


% write source file
fid = fopen([mesh_fn,'.source'],'w');
for i = 1:16
    fprintf(fid,'%f %f \n',sd(i,1),sd(i,2));
end
fclose(fid);

% write meas file
fid = fopen([mesh_fn,'.meas'],'w');
for i = 1:16
    fprintf(fid,'%f %f \n',sd(i,1),sd(i,2));
end
fclose(fid);
% fprintf(fid,'%f %f %f \n',sd(i,1),sd(i,2),sd(i,3));
close all;
clear fid a cor_bmp cor_mim dist f fid_pts i j mesh move1 move2 rho1 rho2
clear sd sd_fn shift str theta1 theta2 tline trans x y z ans
mesh = load_mesh(mesh_fn);
plotmesh(mesh,1)