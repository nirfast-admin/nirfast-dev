% Goes through entire mesh creation process from mimics.

mimics_fn = 'homog_11k';
sd_fn = 'homog_cw_sd.txt';
mesh_fn = 'homog_cw';

mesh = mimics2nirfast_3dmesh_mike(mimics_fn,'tempmesh');
l = importdata('link_full.txt');
fid = fopen('tempmesh.link','w');
for i = 1:16
    fprintf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',...
        l(i,1),l(i,2),l(i,3),l(i,4),l(i,5),l(i,6),l(i,7),l(i,8),...
        l(i,9),l(i,10),l(i,11),l(i,12),l(i,13),l(i,14),l(i,15));
end
fclose(fid)
disp('Placing sources and detectors...') 
mesh = place_sd(sd_fn,'tempmesh');
%%
mesh = minband_opt(mesh);
disp('Moving sources and detectors...')
mesh = move_sd_3d_slab(mesh,4,'opposingfibers.txt');
mesh.excoef = [mesh.excoef(1:2,:);mesh.excoef(4:7,:)];
mesh.wv = [mesh.wv(1:2);mesh.wv(4:7)];
save_mesh(mesh,mesh_fn);

clear ans fid i l mesh_fn mimics_fn sd_fn