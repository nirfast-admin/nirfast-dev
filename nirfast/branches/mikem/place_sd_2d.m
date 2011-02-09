% Places sources and detectors for spec/tomo system.
% Must have SD txt file to load with 4 points at the locations of the
% fiducial markers from the fiber holder.
% M. Mastanduno 6/20/10


% load fiducial points
sd_fn = 'plane2_sd.txt';
mesh_fn = 'temp'; 
pix_sp = .664;

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
% account for bitmap transformation
y = (512-(abs(y)/pix_sp))*pix_sp;
f = [x y z];

dist = 8.9 + (0:12.7:7*12.7);

% place 1-8
shift(1:8,:) = repmat(f(1,1:3),8,1);
move1 = f(2,1:3) - f(1,1:3);
[sd(1:8,1),r,sd(1:8,3)] = cart2pol(move1(1),move1(2),move1(3));
sd(1:8,2) = dist;
% place 9-16
shift(9:16,:) = repmat(f(4,1:3),8,1);
move2 = f(3,1:3) - f(4,1:3);
[sd(9:16,1),r,sd(9:16,3)] = cart2pol(move2(1),move2(2),move2(3));
%th916 = atand(shift(9,2)/shift(9,1));
%sd(9:16,1) = th916;
sd(16:-1:9,2) = dist;
[sd(:,1) sd(:,2),sd(:,3)] = pol2cart(sd(:,1),sd(:,2),sd(:,3));
sd = sd+shift;

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
close all
mesh = load_mesh(mesh_fn);
plotmesh(mesh,1)