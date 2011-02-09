function mesh = MIMICS_to_NIRFAST_3d_MRI_mouse(fn_mimics_mesh, fn_mimics_sd, fn_nirfast_mesh, w)


% Creates 3-D NIRFAST mesh from files exported from Mimics.  

% 1)  .node, .element and .region files are generated from the ANSYS
% files exported from Mimics.  

% 2)  A generic .param file is created.

% 3)  .source and .meas files are also created from exported MedCAD points (Mimics).  

%     Before saving .source and .meas files, 
%     source-detector positions are moved in the following way:
% 		a) Fit plane to the s-d coordinates saved in fn_mimics_sd
% 		b) Project s-d coordinates into this fitted plane
%       c) Move opposing fibers towards one another.  Sources are
%           moved 1 scattering distance inside the tissue surface, 
%           and detectors are moved just inside the mesh.  These positions 
%           are saved in .source and .meas files with 'Fixed' positions 
%           such that load_mesh will not move them automatically (move_source 
%           and move_detector in load_mesh will position s-d incorrectly).
%
%*****************************************************************
% INPUTS:
% 	fn_mimics_mesh = the root name of the export ANSYS files.  
%       For example, if your ANSYS files are 'mouse1.inp' and 'mouse1.txt',
%       fn_mimics_mesh = 'mouse1'
% 	fn_mimics_sd = full file name of the text file (including .txt)
%       containing the MedCAD points which represent s-d locations (in
%       order, starting with source 1).  
% 	fn_nirfast_mesh = the file name of the mesh you want to save.
% 	w = width of Gaussian source (default = 3mm).
% 
% S Davis 2/2009

%****************************************************************


% default w
if nargin == 3
    w = 3;
end

% Create .node, .element, .param, .region files
mimics2nirfast_3dmesh(fn_mimics_mesh, fn_nirfast_mesh);

% Create .source and .meas files (saves .source, .meas, and .link file).
% Puts s-d in same plane, but does not move these positions inside the mesh
% for modeling
disp('Loading s-d locations and moving into fit 3-D plane')
mimics2nirfast_sd_3d_mouse(fn_mimics_sd, fn_nirfast_mesh);

disp('Loading mesh files and moving sources and detectors within mesh boundary')
% Adjust positions of sources and detectors
mesh = move_sd_3d_mouse(fn_nirfast_mesh, w);

disp('Saving mesh')
% Save final mesh
save_mesh(mesh, fn_nirfast_mesh)