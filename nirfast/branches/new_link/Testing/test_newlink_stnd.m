function err = test_newlink_stnd(mesh_stnd)
%***********************************************************
% Standard
% Check that re-ordering source.coord does not affect result
err = [];

mesh = load_mesh(mesh_stnd);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh_anom = add_blob(mesh,blob);

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);

% check FD stnd
data = femdata(mesh_anom,100);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh,[30 30],100,data,2,lambda,'ordered_all_stnd',0);
A = read_solution(mesh,'ordered_all_stnd');

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh2,[30 30],100,data,2,lambda,'unordered_det_IGmesh_stnd',0);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check FD region stnd
mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_region(mesh,100,data,2,lambda,'ordered_all_stnd',0,[0 1]);
A = read_solution(mesh,'ordered_all_stnd');

mesh2.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_region(mesh2,100,data,2,lambda,'unordered_det_IGmesh_stnd',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check CW stnd
data = femdata(mesh_anom,0);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh,[30 30],data,2,lambda,'ordered_all_stnd',0);
A = read_solution(mesh,'ordered_all_stnd');

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh2,[30 30],data,2,lambda,'unordered_det_IGmesh_stnd',0);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check CW region stnd
mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh,data,2,lambda,'ordered_all_stnd',0,[0 1]);
A = read_solution(mesh,'ordered_all_stnd');

mesh2.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh2,data,2,lambda,'unordered_det_IGmesh_stnd',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_stnd);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh_anom = add_blob(mesh,blob);

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);

% check FD stnd
data = femdata(mesh_anom,100);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh,[30 30],100,data,2,lambda,'ordered_all_stnd',0);
A = read_solution(mesh,'ordered_all_stnd');

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh2,[30 30],100,data,2,lambda,'unordered_det_IGmesh_stnd',0);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check FD region stnd
mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_region(mesh,100,data,2,lambda,'ordered_all_stnd',0,[0 1]);
A = read_solution(mesh,'ordered_all_stnd');

mesh2.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_region(mesh2,100,data,2,lambda,'unordered_det_IGmesh_stnd',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check CW stnd
data = femdata(mesh_anom,0);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh,[30 30],data,2,lambda,'ordered_all_stnd',0);
A = read_solution(mesh,'ordered_all_stnd');

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh2,[30 30],data,2,lambda,'unordered_det_IGmesh_stnd',0);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% check CW region stnd
mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh,data,2,lambda,'ordered_all_stnd',0,[0 1]);
A = read_solution(mesh,'ordered_all_stnd');

mesh2.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh2,data,2,lambda,'unordered_det_IGmesh_stnd',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_stnd);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,100);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data1.link=data.link(dat,:);
data1.paa=data.paa(dat,:);

% check FD stnd and FD region stnd
lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh,[30 30],100,data1,2,lambda,'odds_only_stnd',0);
A = read_solution(mesh,'odds_only_stnd');

mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_region(mesh,100,data1,2,lambda,'odds_only_hp_stnd',0,[0 1]);
Ahp = read_solution(mesh,'odds_only_hp_stnd');

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
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh2_anom = add_blob(mesh2,blob);

data2 = femdata(mesh2_anom,100);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd(mesh2,[30 30],100,data2,2,lambda,'odds_only_2_stnd',0);
B = read_solution(mesh,'odds_only_2_stnd');

mesh2.region = mesh2_anom.region;
lambda=1;
reconstruct_stnd_region(mesh2,100,data2,2,lambda,'odds_only_2hp_stnd',0,[0 1]);
Bhp = read_solution(mesh,'odds_only_2hp_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = [Bhp.mua - Ahp.mua, Bhp.mus - Ahp.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

%  check CW stnd and CW region stnd
mesh = load_mesh(mesh_stnd);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,0);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data1.link=data.link(dat,:);
data1.paa=data.paa(dat,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh,[30 30],data1,2,lambda,'odds_only_stnd',0);
A = read_solution(mesh,'odds_only_stnd');

mesh.region = mesh_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh,data1,2,lambda,'odds_only_hp_stnd',0,[0 1]);
Ahp = read_solution(mesh,'odds_only_hp_stnd');

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
blob.region=1;
blob.mua=0.02;
blob.mus=1.4;
mesh2_anom = add_blob(mesh2,blob);

data2 = femdata(mesh2_anom,0);

lambda.type='Automatic';
lambda.value=10;
reconstruct_stnd_cw(mesh2,[30 30],data2,2,lambda,'odds_only_2_stnd',0);
B = read_solution(mesh,'odds_only_2_stnd');

mesh2.region = mesh2_anom.region;
lambda=1;
reconstruct_stnd_cw_region(mesh2,data2,2,lambda,'odds_only_2hp_stnd',0,[0 1]);
Bhp = read_solution(mesh,'odds_only_2hp_stnd');

diff = [B.mua - A.mua, B.mus - A.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

diff = [Bhp.mua - Ahp.mua, Bhp.mus - Ahp.mus];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

