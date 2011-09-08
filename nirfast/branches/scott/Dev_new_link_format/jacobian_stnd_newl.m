function [J,data,mesh]=jacobian_stnd_newl(mesh,frequency,mesh2)

% [J,data,mesh]=jacobian_stnd(mesh,frequency,mesh2)
%
% Calculates the Jacobian (both complex version and separate parts
% in terms of kappa and mua for log amplitude and phase
% (radians). See any of Dartmouth Publications regarding the
% structure. Also calculates data (phase and amplitude)
% outputs phase and amplitude in structure data
% and mesh information in mesh
% 
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% mesh2 is optional mesh basis

source = unique(mesh.link(:,1));

if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

% modulation frequency
omega = 2*pi*frequency*1e6;

% Create FEM matricex
if mesh.dimension == 2
  [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
			    sort(mesh.elements')', ...
			    mesh.bndvtx,...
			    mesh.mua,...
			    mesh.kappa,...
			    mesh.ksi,...
			    mesh.c,...
			    omega);
elseif mesh.dimension ==3
  [i,j,s] = gen_matrices_3d(mesh.nodes,...
			    sort(mesh.elements')', ...
			    mesh.bndvtx,...
			    mesh.mua,...
			    mesh.kappa,...
			    mesh.ksi,...
			    mesh.c,...
			    omega);
end

junk = length(find(i==0));
MASS = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
clear junk i j s omega

% If the fn.ident exists, then we must modify the FEM matrices to
% account for refractive index mismatch within internal boundaries
if isfield(mesh,'ident') == 1
  disp('Modifying for refractive index')
  M = bound_int(MASS,mesh);
  MASS = M;
  clear M
end

% Calculate the RHS (the source vectors. For simplicity, we are
% just going to use a Gaussian Source, The width of the Gaussian is
% changeable (last argument). The source is assumed to have a
% complex amplitude of complex(cos(0.15),sin(0.15));

% Now calculate source vector
% NOTE last term in mex file 'qvec' is the source FWHM
%
[nnodes,junk]=size(mesh.nodes);
[nsource,junk]=size(source);
qvec = spalloc(nnodes,nsource,nsource*100);
if mesh.dimension == 2
  for i = 1 : nsource
    if mesh.source.fwhm(i) == 0
        qvec(:,i) = gen_source_point(mesh,mesh.source.coord(source(i),1:2));
    else
      qvec(:,i) = gen_source(mesh.nodes(:,1:2),...
			   sort(mesh.elements')',...
			   mesh.dimension,...
			   mesh.source.coord(source(i),1:2),...
			   mesh.source.fwhm(source(i)));
    end
  end
elseif mesh.dimension == 3
  for i = 1 : numel(source)
    if mesh.source.fwhm(i) == 0
        qvec(:,i) = gen_source_point(mesh,mesh.source.coord(source(i),1:3));
    else
    qvec(:,i) = gen_source(mesh.nodes,...
			   sort(mesh.elements')',...
			   mesh.dimension,...
			   mesh.source.coord(source(i),:),...
			   mesh.source.fwhm(source(i)));
    end
  end
end
clear junk i nnodes nsource w;

% Catch zero frequency (CW) here
if frequency == 0
  MASS = real(MASS);
  qvec = real(qvec);
end

% catch error in source vector
junk = sum(qvec);
junk = find(junk==0);
if ~isempty(junk)
    display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
end
clear junk

% Calculate field for all sources
[data.phi,mesh.R]=get_field(MASS,mesh,qvec);
clear qvec;

% Now calculate Adjoint source vector
[qvec] = gen_source_adjoint_newl(mesh);

% Catch zero frequency (CW) here
if frequency == 0
  qvec = real(qvec);
end

% Calculate adjoint field for all detectors
[data.aphi]=get_field(conj(MASS),mesh,conj(qvec));
clear qvec MASS;

% Calculate boundary data
[data.complex]=get_boundary_data_newl(mesh,data.phi);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);

data.phase = atan2(imag(data.complex),...
		   real(data.complex));
data.phase(find(data.phase<0)) = data.phase(find(data.phase<0)) + (2*pi);
data.phase = data.phase*180/pi;

data.paa = [data.amplitude data.phase];

%data.phi(find(isnan(data.phi)==1))=[];
%data.aphi(find(isnan(data.aphi)==1))=[];

if nargin == 3 % use second mesh basis for jacobian
    data2 = interpolatef2r(mesh,mesh2,data);
    data2.complex = data.complex;
    ind = find(data.link(:,3) == 0);
    data2.complex(ind,:)=[];
    
    % Calculate Jacobian
    % Catch zero frequency (CW) here
    if frequency == 0
        [J] = build_jacobian_cw(mesh2,data2);
    else
        [J] = build_jacobian_newl(mesh2,data2);
    end
elseif nargin == 2
    data2.complex = data.complex;
    ind = find(data.link(:,3) == 0);
    data2.complex(ind,:)=[];
    
    % Calculate Jacobian
    % Catch zero frequency (CW) here
    if frequency == 0
        [J] = build_jacobian_cw(mesh,data2);
    else
        [J] = build_jacobian_newl(mesh,data2);
    end
end


function data2 = interpolatef2r(fwd_mesh,recon_mesh,data)
% This function interpolates fwd_mesh data into recon_mesh
% Used to calculate the Jacobian on second mesh

for i = 1 : length(recon_mesh.nodes)
    if fwd_mesh.fine2coarse(i,1) ~= 0
    data2.phi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
    data.phi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
    data2.aphi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
    data.aphi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
    elseif fwd_mesh.fine2coarse(i,1) == 0
    dist = distance(mesh.nodes,...
                    mesh.bndvtx,...
                    pixel.nodes(i,:));
    mindist = find(dist==min(dist));
    mindist = mindist(1);
    data2.phi(i,:) = data.phi(mindist,:);
    data2.aphi(i,:) = data.phi(mindist,:);
    end
end