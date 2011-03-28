function err = test_new_link_spec(mesh_spec)
err = [];
% spectral
% Check that re-ordering source.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,100);
save_data(data,'data_ordered_sd.paa');
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],100,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],100,data,2,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,100);
save_data(data,'data_ordered_sd.paa');
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],100,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],100,data,2,lambda,'unordered_det_IGmesh',0);
B = read_solution(mesh,'unordered_det_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,100);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],100,data,2,lambda,'odds_only',0);
A = read_solution(mesh,'odds_only');

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
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh2_anom = add_blob(mesh2,blob);
save_mesh(mesh2_anom,'mesh_anom');
data = femdata(mesh2_anom,100);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],100,data,2,lambda,'odds_only_2',0);
B = read_solution(mesh,'odds_only_2');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end


% CW spectral*****************************************************
% Check that re-ordering source.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,0);
save_data(data,'data_ordered_sd.paa');
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,0);
save_data(data,'data_ordered_sd.paa');
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);
lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'unordered_det_IGmesh',0);
B = read_solution(mesh,'unordered_det_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
save_mesh(mesh_anom,'mesh_anom');
data = femdata(mesh_anom,0);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh,[30 30],0,data,2,lambda,'odds_only',0);
A = read_solution(mesh,'odds_only');

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
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh2_anom = add_blob(mesh2,blob);
save_mesh(mesh2_anom,'mesh_anom');
data = femdata(mesh2_anom,0);

lambda.type='Automatic';
lambda.value=10;
reconstruct(mesh2,[30 30],0,data,2,lambda,'odds_only_2',0);
B = read_solution(mesh,'odds_only_2');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end


% Spectral region*****************************************************
% Check that re-ordering source.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,100);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_region(mesh,100,data,2,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_region(mesh2,100,data,2,lambda,'unordered_source_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
data = femdata(mesh_anom,100);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_region(mesh,100,data,2,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_region(mesh2,100,data,2,lambda,'unordered_det_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
data = femdata(mesh_anom,100);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_region(mesh,100,data,2,lambda,'odds_only',0,[0 1]);
A = read_solution(mesh,'odds_only');

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
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh2_anom = add_blob(mesh2,blob);
data = femdata(mesh2_anom,100);

mesh2.region = mesh2_anom.region;
lambda=10;
reconstruct_spectral_region(mesh2,100,data,2,lambda,'odds_only_2',0,[0 1]);
B = read_solution(mesh,'odds_only_2');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Spectral CW region*****************************************************
% Check that re-ordering source.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);

data = femdata(mesh_anom,0);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh,data,2,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.source.coord(1:3,:) = mesh.source.coord(5:7,:);
mesh2.source.coord(5:7,:) = mesh.source.coord(1:3,:);
mesh2.source.num(1:3,:) = mesh.source.num(5:7,:);
mesh2.source.num(5:7,:) = mesh.source.num(1:3,:);

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh2,data,2,lambda,'unordered_source_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that re-ordering meas.coord does not affect result
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
data = femdata(mesh_anom,0);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh,data,2,lambda,'ordered_all',0,[0 1]);
A = read_solution(mesh,'ordered_all');

mesh2 = mesh;
mesh2.meas.coord(1:3,:) = mesh.meas.coord(5:7,:);
mesh2.meas.coord(5:7,:) = mesh.meas.coord(1:3,:);
mesh2.meas.num(1:3,:) = mesh.meas.num(5:7,:);
mesh2.meas.num(5:7,:) = mesh.meas.num(1:3,:);
mesh2.meas.int_func(1:3,:) = mesh.meas.int_func(5:7,:);
mesh2.meas.int_func(5:7,:) = mesh.meas.int_func(1:3,:);

mesh2.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh2,data,2,lambda,'unordered_det_IGmesh',0,[0 1]);
B = read_solution(mesh,'unordered_det_IGmesh');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end

% Check that using only odd numbers in data file without place holders
% matches leaving these out from the begininng
mesh = load_mesh(mesh_spec);
blob.x=20;
blob.y=10;
blob.r=10;
blob.z=0;
blob.region=1;
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh_anom = add_blob(mesh,blob);
data = femdata(mesh_anom,0);

dat = logical(mod(data.link(:,1),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);
dat = logical(mod(data.link(:,2),2));
data.link=data.link(dat,:);
data.paa=data.paa(dat,:);

mesh.region = mesh_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh,data,2,lambda,'odds_only',0,[0 1]);
A = read_solution(mesh,'odds_only');

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
blob.sa=1.2;
blob.sp=1.3;
blob.HbO=0.02;
blob.deoxyHb=0.015;
blob.Water=0.8;
mesh2_anom = add_blob(mesh2,blob);
data = femdata(mesh2_anom,0);

mesh2.region = mesh2_anom.region;
lambda=10;
reconstruct_spectral_cw_region(mesh2,data,2,lambda,'odds_only_2',0,[0 1]);
B = read_solution(mesh,'odds_only_2');

diff = [B.conc - A.conc, B.sp - A.sp, B.sa - A.sa];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end
