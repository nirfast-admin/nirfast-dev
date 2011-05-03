function testing
% Set optical properties in mesh.
fn = 'sdh4_session1_3region';
fn1 = 'SDH4';

%[mesh] = load_mesh(['uncal_',fn]);
%mesh = set_mua_mus_mousebrain_09292009(mesh);
%save_mesh(mesh,['uncal_',fn])
%delete(['uncal_',fn,'.link']);


% AF750
mesh = load_mesh('SDH4_session1_3region');
mesh = set_mua_mus_mousebrain_09292009(mesh);
frame = 1;
fndata = ['CAL_',fn1,'_immediate_post_AF750_frame'];
%fndata_save = ['CAL4recon_',fn,'_AF750_frame'];
%fnIG_save = ['IG_',fn,'_AF750_frame'];
for i = 1:length(frame)
    tic
    [data_cal,IG_mesh]=calibrate_fl_secant2(mesh,[fndata,num2str(i),'.paa'],0,50,10e-4);
    toc
    tic
    [data_cal,IG_mesh]=calibrate_fl(mesh,[fndata,num2str(i),'.paa'],0,50,10e-4);
    toc
    %save_data([fndata_save,num2str(i),'.paa'],data_cal);
    %save_mesh(IG_mesh,[fnIG_save,num2str(i)]);
end