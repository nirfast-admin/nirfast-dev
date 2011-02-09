function mesh = place_sd(sd_fn,mesh_fn)

% Places sources and detectors for spec/tomo system.
% Must have SD txt file to load with 4 points at the locations of the
% fiducial markers from the fiber holder.
% M. Mastanduno 6/20/10


% load fiducial points
% sd_fn = 'ph92110-3d_p3.txt';
% mesh_fn = 'p3'; 

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
f = [x y z];

dist = 8.9 + (0:12.7:7*12.7);

% place 1-8
shift(1:8,:) = repmat(f(1,1:3),8,1);
move1 = f(2,1:3) - f(1,1:3);
if abs(y(1)-y(2))<.1
    [sd(1:8,1),r,sd(1:8,3)] = cart2pol(move1(1),move1(3),move1(2));
else
    [sd(1:8,1),r,sd(1:8,3)] = cart2pol(move1(1),move1(2),move1(3));
end
sd(1:8,2) = dist;
% place 9-16
shift(9:16,:) = repmat(f(4,1:3),8,1);
move2 = f(3,1:3) - f(4,1:3);
if abs(y(1)-y(2))<.1
    [sd(9:16,1),r,sd(9:16,3)] = cart2pol(move2(1),move2(3),move2(2));
else
    [sd(9:16,1),r,sd(9:16,3)] = cart2pol(move2(1),move2(2),move2(3));
end
sd(16:-1:9,2) = dist;
if abs(y(1)-y(2))<.1
    [sd(:,1) sd(:,3),sd(:,2)] = pol2cart(sd(:,1),sd(:,2),sd(:,3));
else
    [sd(:,1) sd(:,2),sd(:,3)] = pol2cart(sd(:,1),sd(:,2),sd(:,3));
end
sd = sd+shift;


%sdx = [sdx18 sdx916]';
%sdy = [sdz18 sdz916]';
%sd = [sdx sdy sdz];

% write source file
fid = fopen([mesh_fn,'.source'],'w');
fprintf(fid,'fixed\n');
for i = 1:16
    fprintf(fid,'%f %f %f \n',sd(i,1),sd(i,2),sd(i,3));
end
fclose(fid);

% write meas file
fid = fopen([mesh_fn,'.meas'],'w');
fprintf(fid,'fixed\n');
for i = 1:16
    fprintf(fid,'%f %f %f \n',sd(i,1),sd(i,2),sd(i,3));
end
fclose(fid);

close all
mesh = load_mesh(mesh_fn);
plotmesh(mesh,1)
end