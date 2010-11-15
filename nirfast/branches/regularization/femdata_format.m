function ref = femdata_format(mesh,kappa,mua,frequency)
% mesh is the geometry info
% kappa is the vector of diffusion
% mua is the vector of absorption
% frequency is the modulation frequency
% ref is the resulting data, with amplitude on the top
%   half and phase on the bottom half

mesh.kappa = kappa;
mesh.mua = mua;
mesh.mus = (1./(3.*mesh.kappa))-mesh.mua;
ref = femdata(mesh,frequency);
ref(:,1) = log(data.amplitude);
ref(:,2) = data.phase;
ref(:,2) = ref(:,2)/180.0*pi;
ref(find(ref(:,2)<0),2) = ref(find(ref(:,2)<0),2) + (2*pi);
ref(find(ref(:,2)>(2*pi)),2) = ref(find(ref(:,2)>(2*pi)),2) - (2*pi);
ref = reshape(ref',length(ref)*2,1);