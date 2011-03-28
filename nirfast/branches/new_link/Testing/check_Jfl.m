function check_Jfl(mesh_fl)
% Test whether normalization of Jm is necessary.
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


lambda.type='Automatic';
lambda.value=10;
reconstruct_fl(mesh,[30 30],0,data,10,lambda,'ordered_all',0);
A = read_solution(mesh,'ordered_all');

lambda.type='Automatic';
lambda.value=10;
reconstruct_fl2(mesh,[30 30],0,data,10,lambda,'unordered_source_IGmesh',0);
B = read_solution(mesh,'unordered_source_IGmesh');

diff = [B.etamuaf - A.etamuaf];
diff = sum(sum(diff));
if diff == 0
    disp('OK!')
    err = [err; 'OK!'];
else
    err = [err; 'E>0'];
end