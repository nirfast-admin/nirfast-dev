function test_femdata_fl

mesh = load_mesh('circle2000_86_fl');

tic
% calc mm, fl, x
data1 = femdata(mesh,0)
toc

tic
% calc fl, x
mesh.mm = 0;
data2 = femdata(mesh,0)
toc

tic
% calc mm, fl, x
mesh.mm = 1;
data3 = femdata(mesh,0)
toc

tic
% calc mm, fl, x
mesh.fl = 1;
data4 = femdata(mesh,0)
toc

tic
% calc mm, x
mesh.fl = 0;
data5 = femdata(mesh,0)
toc

tic
% calc mm
mesh.phix = data1.phix;
data6 = femdata(mesh,0)
toc

tic
% calc mm, fl
mesh.fl = 1;
data7 = femdata(mesh,0)
toc

tic
% calc fl
mesh.mm = 0;
data8 = femdata(mesh,0)
toc

% Calibrate and Recon
circle2000_86_fl = load_mesh('C:\NewNirfast\meshes\fl\circle2000_86_fl');
blob.x=0;
blob.y=20;
blob.r=8;
blob.z=0;
blob.muaf=0.005;
blob.region=1;
circle2000_86_fl_anom = add_blob(circle2000_86_fl,blob);
circle2000_86_fl_anom_data = femdata(circle2000_86_fl_anom,0);
% New femdata
tic
[data_cal,mesh_cal]=calibrate_fl(circle2000_86_fl,circle2000_86_fl_anom_data,0,50,10e-5);
toc

% Old version
tic
[data_cal,mesh_cal]=calibrate_fl_old(circle2000_86_fl,circle2000_86_fl_anom_data,0,50,10e-5);
toc

% Result:  calibrate_fl runs in 45% the time of the old version - a 55%
% improvement!

% Checking Jacobian change at each iteration
lambda.type='Automatic';
lambda.value=10;
tic
[mesh,pj] = reconstruct_fl(circle2000_86_fl,[30 30],0,circle2000_86_fl_anom_data,40,lambda,'junk_test1',0);
toc
read_solution(mesh,'junk_test1');

tic
[mesh,pj] = reconstruct_fl2(circle2000_86_fl,[30 30],0,circle2000_86_fl_anom_data,40,lambda,'junk_test2',0);
toc
read_solution(mesh,'junk_test2');