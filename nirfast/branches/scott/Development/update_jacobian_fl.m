function [Jout,data]=update_jacobian_fl(mesh,MASS_m,frequency,datax,Jin)

%  Comments


% error checking
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

% Calculate the RHS (the source vectors) for the Emission.
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);
[nsource,junk]=size(source);
qvec = zeros(nnodes,nsource);
% Simplify the RHS of emission equation
beta = mesh.gamma.*(1-(sqrt(-1).*omega.*mesh.tau));
% get rid of any zeros!
if frequency == 0
    beta(beta==0) = 1e-20;
else
    beta(beta==0) = complex(1e-20,1e-20);
end

if mesh.dimension == 2
    for i = 1 : nsource
        val = beta.*datax.phi(:,i);
        qvec(:,i) = gen_source_fl(mesh.nodes(:,1:2),...
            sort(mesh.elements')',...
            mesh.dimension,...
            val);
    end
elseif mesh.dimension == 3
    for i = 1 : nsource
        val = beta.*datax.phi(:,i);
        qvec(:,i) = gen_source_fl(mesh.nodes,...
            sort(mesh.elements')',...
            mesh.dimension,...
            val);
    end
end

clear junk i nnodes nsource val beta;

% Calculate EMISSION field for all sources
[data.phim,mesh.R]=get_field(MASS_m,mesh,qvec);
clear qvec;

% Calculate boundary data
[data.complexm]=get_boundary_data(mesh,data.phim);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitudem = abs(data.complexm);

data.phasem = atan2(imag(data.complexm),...
    real(data.complexm));
data.phasem(data.phasem<0) = data.phasem(data.phasem<0) + (2*pi);
data.phasem = data.phasem*180/pi;

data.paam = [data.amplitudem data.phasem];
data.phix = datax.phi;

% Update the Emission jacobian
data2 = data;
ind = data.link(:,3) == 0;
data2.complexm(ind,:)=[];
[nsd, msd] = size(mesh.link);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        Jout.completem(k,:) = ...
            real(Jin.complexm(k,1:end/2)./data2.complexm(k));
    end
    k = k + 1;
end