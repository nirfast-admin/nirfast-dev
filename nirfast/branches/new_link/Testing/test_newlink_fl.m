function err = test_newlink_fl(mesh_fl)

%***********************************************************
% Fluorescence
% Check that re-ordering source.coord does not affect result
err = [];

mesh = load_mesh(mesh_fl);

blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.muaf=0.005;
blob.region=1;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,0);

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_fl_region(mesh,0,data,3,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_fl_region(mesh2,0,data,3,lambda,'unordered_source_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh,[30 30; 10 10],0,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh2,[30 30; 10 10],0,data,2,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

[datajunk1,meshjunk1] = calibrate_fl(mesh, data,...
    0, 50, 10^-3)

[datajunk2,meshjunk2] = calibrate_fl(mesh2, data,...
    0, 50, 10^-3)

diff = meshjunk1.muaf - meshjunk2.muaf;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = datajunk1.amplitudefl - datajunk2.amplitudefl;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

%*******
% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_fl);

blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.muaf=0.005;
blob.region=1;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,0);

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'ordered_all_fl',0);
A = read_solution(mesh,'ordered_all_fl');

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'unordered_det_IGmesh_fl',0);
B = read_solution(mesh,'unordered_det_IGmesh_fl');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end


mesh.region = mesh_anom.region;
lambda=10;
reconstruct_fl_region(mesh,0,data,3,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_fl_region(mesh2,0,data,3,lambda,'unordered_source_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh,[30 30; 10 10],0,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh2,[30 30; 10 10],0,data,2,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

[datajunk1,meshjunk1] = calibrate_fl(mesh, data,...
    0, 50, 10^-3)

[datajunk2,meshjunk2] = calibrate_fl(mesh2, data,...
    0, 50, 10^-3)

diff = meshjunk1.muaf - meshjunk2.muaf;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = datajunk1.amplitudefl - datajunk2.amplitudefl;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end
%*******
% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_fl);

blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.muaf=0.005;
blob.region=1;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,0);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paafl=data.paafl(dat,:);
data.amplitudefl=data.amplitudefl(dat,:);
data.amplitudex=data.amplitudex(dat,:);
dat = logical(mod(data.link(:,2),2));
data.link=data.link(dat,:);
data.paafl=data.paafl(dat,:);
data.amplitudefl=data.amplitudefl(dat,:);
data.amplitudex=data.amplitudex(dat,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'odds_only_fl',0);
A = read_solution(mesh,'odds_only_fl');

mesh.region = mesh_anom.region;
lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh,[30 30],0,data,2,lambda,'odds_only_fl',0);
Aspat = read_solution(mesh,'odds_only_fl');

lambda=10;
reconstruct_fl_region(mesh,0,data,3,lambda,'ordered_all',0,[0 1]);
Ahp = read_solution(mesh,'ordered_all');

[datajunk1,meshjunk1] = calibrate_fl(mesh, data,...
    0, 50, 10^-3)

mesh2 = mesh;
mesh2.source.num(2:2:end,:)=[];
mesh2.source.coord(2:2:end,:)=[];
mesh2.source.fwhm(2:2:end,:)=[];
mesh2.meas.num(2:2:end,:)=[];
mesh2.meas.int_func(2:2:end,:)=[];
mesh2.meas.coord(2:2:end,:)=[];

dat = logical(mod(mesh2.link(:,1),2));
mesh2.link=mesh2.link(dat,:);
dat = logical(mod(mesh2.link(:,2),2));
mesh2.link=mesh2.link(dat,:);

blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.muaf=0.005;
blob.region=1;
mesh2_anom = add_blob(mesh2,blob);

data = femdata(mesh2_anom,0);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'odds_only_fl_2',0);
B = read_solution(mesh,'odds_only_fl_2');

mesh2.region = mesh_anom.region;
lambda.type='Automatic';
lambda.value=10;
reconstruct_fl_spatial(mesh2,[30 30],0,data,2,lambda,'odds_only_fl_2',0);
Bspat = read_solution(mesh,'odds_only_fl_2');

lambda=10;
reconstruct_fl_region(mesh2,0,data,3,lambda,'unordered_source_IGmesh',0,[0 1]);
Bhp = read_solution(mesh,'unordered_source_IGmesh');

[datajunk2,meshjunk2] = calibrate_fl(mesh2, data,...
    0, 50, 10^-3)

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = [Bspat.etamuaf - Aspat.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = [Bhp.etamuaf - Ahp.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = meshjunk1.muaf - meshjunk2.muaf;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = datajunk1.amplitudefl - datajunk2.amplitudefl;
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end
