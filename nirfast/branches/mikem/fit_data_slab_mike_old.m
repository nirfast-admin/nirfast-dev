function [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data_slab_mike(mesh,...
                                                               data,...
                                                               frequency,...
                                                               iteration)
                                                           
% [mua,mus,lnI_offset,phase_offset,fem_data] = fit_data(mesh,...
%                                   data,frequency,iteration)
%
% fits data to a given model to give initial guesses for
% reconstruction as well as data offsets
%
% mesh is the input mesh (variable)
% data is the boundary data (variable)
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% outputs are the intial guesses

% load link index from a text file. 0's are for reflectance, 1's for
% transmission
link_key = importdata('link_reflectance_transmission.txt','\t',3);
link_key = link_key.data;
% calculate the source / detector distance for each measurement
k = 1; r = 1; t = 1;
datanum = 0;
[ns,junk]=size(mesh.source.coord);
for i = 1 : ns
  for j = 1 : length(mesh.link(i,:))
      datanum = datanum + 1;
      % deal with nans here
      if isnan(data(datanum,1)) || isnan(data(datanum,2))
          mesh.link(i,j) = 0;
      end
      % if link is nonzero
      if mesh.link(i,j) ~= 0
          jj = mesh.link(i,j);
          if link_key(i,j) == 0 % reflectance only
              dist_r(r,1) = sqrt(sum((mesh.source.coord(i,:) - ...
                  mesh.meas.coord(jj,:)).^2));
              data_r(r,:) = data(datanum,:);
              index_meas(k) = 0;
              r = r + 1;
          elseif link_key(i,j) == 1 % transmission only
              dist_t(t,1) = sqrt(sum((mesh.source.coord(i,:) - ...
                  mesh.meas.coord(jj,:)).^2));
              data_t(t,:) = data(datanum,:);
              index_meas(k) = 1;
              t = t + 1;
          end
          % still get all distances
          dist(k,1) = sqrt(sum((mesh.source.coord(i,:) - ...
              mesh.meas.coord(jj,:)).^2));
          k = k+1;
      else
        data(datanum,1) = NaN;
      end
  end
end

% % convert log amplitude into amplitude
% data(:,1) = exp(data(:,1));

% Set lnrI, lnr and phase!
[j,k] = size(data(:,1));
[j2,k2] = size(dist);

% deal with NaNs in all 3 modes separately. first total
dist_orig = dist;
ind = unique([find(isnan(data(:,1))==1); find(isnan(data(:,2))==1)]);
ind = setdiff(1:size(data,1),ind);
data = data(ind,:);
% now reflectance
dist_orig_r = dist_r;
ind_r = unique([find(isnan(data_r(:,1))==1); find(isnan(data_r(:,2))==1)]);
ind_r = setdiff(1:size(data_r,1),ind_r);
data_r = data_r(ind_r,:);
% now transmission
dist_orig_t = dist_t;
ind_t = unique([find(isnan(data_t(:,1))==1); find(isnan(data_t(:,2))==1)]);
ind_t = setdiff(1:size(data_t,1),ind_t);
data_t = data_t(ind_t,:);

lnrI_ref = log(data_r(:,1).*dist_r);
lnrI_trans = log(data_t(:,1).*dist_t);
lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));

phase = data(:,2);
phase_r = data_r(:,2);

% make plots for ref, trans, all
figure;
subplot(1,2,1);
plot(dist_r,lnrI_ref,'.');
ylabel('lnrI Reflectance data');
xlabel('Source / Detector distance');

subplot(1,2,2);
plot(dist_t,lnrI_trans,'.');
ylabel('lnrI Transmittance data');
xlabel('Source / Detector distance');

figure;
subplot(2,2,1);
plot(dist,lnrI,'.')
ylabel('lnrI');
xlabel('Source / Detector distance');
subplot(2,2,2);
plot(dist,phase,'.')
ylabel('Phase');
xlabel('Source / Detector distance');
drawnow
pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase
m0_r = polyfit(dist_r,phase_r,1); m0_r = m0_r(1);
% calculate the slopes for reflectance and transmission seperately
m1_r = polyfit(dist_r,lnrI_ref,1); m1_r = m1_r(1);
disp(['Slope from reflectance only is: ',num2str(m1_r)]);
m1 = m1_r;
m0 = m0_r;
m1_old = polyfit(dist,lnrI,1); m1_old = m1_old(1);

% fit data using an analytical model
% based on Pogue paper
omega = 2*pi*frequency*1e6;
c=(3e11./mean(mesh.ri));

mua = 0.01; mus = 1; kappa = 1/(3*(mua+mus));

for i = 1 : 25
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mua = mua + 0.0001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha1 = -f1*cos(0.5*atan2(omega/c,mua));
  phi1 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mua = mua - 0.0001;   
  da = (alpha-alpha1)/0.0001;
  ds = (phi-phi1)/0.0001;
  mua = mua - (m1-alpha)/da;
  mus = mus - (m0-phi)/ds*0.1;
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha = -f1*cos(0.5*atan2(omega/c,mua));
  phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  
  mus = mus + 0.001; 
  kappa = 1/(3*(mua+mus));
  f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
  alpha2 = -f1*cos(0.5*atan2(omega/c,mua));
  phi2 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
  mus = mus - 0.001; 
  da = (alpha-alpha2)/0.001;
  ds = (phi-phi2)/0.001;
  mua = mua - (m1-alpha)/da*0.01;
  mus = mus - (m0-phi)/ds*0.2;
  kappa = 1/(3*(mua+mus));
end

disp('Global values calculated from Analytical fit');
disp(['Absorption = ' num2str(mua) ' mm-1']);
disp(['Scatter    = ' num2str(mus) ' mm-1']);
disp('============================================');

% Set the global values onto mesh.
mesh.mua(:) = mua;
mesh.mus(:) = mus;
mesh.kappa(:) = kappa;

% Fit for mua and mus using FEM
dist = dist_orig;
jj = 0;
 ind_r = find(index_meas(:) == 0);
 disp('numerical fit ................');
while jj ~= iteration
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
    
  femlnrI_r = log(fem_data(ind_r,1).*dist_r);
  femphase_r = fem_data(ind_r,2);
  
  phi0_r = polyfit(dist_r,femphase_r,1); phi0_r = phi0_r(1);
  alpha0_r = polyfit(dist_r,femlnrI_r,1); alpha0_r = alpha0_r(1);
  
  mesh.mua(:) = mesh.mua(:)+0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI_r = log(fem_data(ind_r,1).*dist_r);
  femphase_r = fem_data(ind_r,2);
  
  phi1_r = polyfit(dist_r,femphase_r,1); phi1_r = phi1_r(1);
  alpha1_r = polyfit(dist_r,femlnrI_r,1); alpha1_r = alpha1_r(1);
  
  mesh.mua(:) = mesh.mua(:)-0.0001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  alpha0 = alpha0_r;
  alpha1 = alpha1_r;
  phi0 = phi0_r;
  phi1 = phi1_r;
  
  da = (alpha0-alpha1)/0.0001;
  ds = (phi0-phi1)/0.0001;
  
  mesh.mua(:) = mesh.mua(:) - (m1-alpha0)/da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.1;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
 [mesh.mua(1) mesh.mus(1)];
 
 %
 
 [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI_r = log(fem_data(ind_r,1).*dist_r);
  femphase_r = fem_data(ind_r,2);

  phi0_r = polyfit(dist_r,femphase_r,1); phi0_r = phi0_r(1);
  alpha0_r = polyfit(dist_r,femlnrI_r,1); alpha0_r = alpha0_r(1);
  alpha0 = alpha0_r;
  phi0 = phi0_r;
  
  mesh.mus(:) = mesh.mus(:)+0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  [fem_data]=femdata(mesh,frequency);
  fem_data = fem_data.paa;
  
  femlnrI_r = log(fem_data(ind_r,1).*dist_r);
  femphase_r = fem_data(ind_r,2);
  
  phi2_r = polyfit(dist_r,femphase_r,1); phi2_r = phi2_r(1);
  alpha2_r = polyfit(dist_r,femlnrI_r,1); alpha2_r = alpha2_r(1);
  alpha2 = alpha2_r;
  phi2 = phi2_r;
  
  mesh.mus(:) = mesh.mus(:)-0.001;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  da = (alpha0-alpha2)/0.001;
  ds = (phi0-phi2)/0.001;
  
  mesh.mua(:) = mesh.mua(:) - abs((m1-alpha0)/da*0.002);
  mesh.mua(find(mesh.mua(:)<0)) = mesh.mua(find(mesh.mua(:)<0)) + da;
  mesh.mus(:) = mesh.mus(:) - (m0-phi0)/ds.*0.2;
  mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
  
  err_a = abs(mean(mesh.mua)-mua)./mua;
  err_s = abs(mean(mesh.mus)-mus)./mus;
  
  if ((err_a < 0.001) & (err_s < 0.001))
    jj = iteration;
  else
    jj = jj + 1;
  end
  
  mua = mean(mesh.mua);
  mus = mean(mesh.mus);
  disp('Global values calculated from Numerical fit');
  disp(['Iteration ' num2str(jj) ' of ' num2str(iteration)]);
  disp(['Absorption = ' num2str(mua) ' mm-1 with error of ' num2str(err_a)]);
  disp(['Scatter    = ' num2str(mus) ' mm-1 with error of ' num2str(err_s)]);
  disp('-------------------------------------------------');
end

% Calculate data based on these global values
[fem_data]=femdata(mesh,frequency);
fem_data = fem_data.paa;


% Arrange data to calculate offset
femlnI = log(fem_data(:,1));
femphase = fem_data(:,2);

% Find offset
% make sure any original NaNs are present
lnI = log(data(:,1));
phase = data(:,2);

% Set offset based on particular source / detector
% we do this because data is not symmetrical!
[n,m]=size(mesh.link);
lnI_offset_tmp = spalloc(n,m,n*m);
phase_offset_tmp = spalloc(n,m,n*m);
k = 1;
for i = 1 : n
  for j = 1 : m
    if mesh.link(i,j) ~= 0
      jj = mesh.link(i,j);
      lnI_offset_tmp(i,jj) = lnI(k)-femlnI(k);
      phase_offset_tmp(i,jj) = phase(k)-femphase(k);
      k = k + 1;
    end
  end
end

%lnI_offset = repmat(mean(lnI_offset')',1,m);
%phase_offset = repmat(mean(phase_offset')',1,m);
lnI_offset = zeros(size(femlnI));
phase_offset = zeros(size(femphase));
k = 1;
for i = 1 : n
  for j = 1 : m
    if mesh.link(i,j) ~= 0
      temp = lnI_offset_tmp(i,:);
      temp(isnan(temp)) = [];
      lnI_offset(k) = mean(temp);
      temp = phase_offset_tmp(i,:);
      temp(isnan(temp)) = [];
      phase_offset(k) = mean(temp);
      k = k + 1;
    end
  end
end

subplot(2,2,3);
plot(lnI,'k');
hold on
plot(femlnI+lnI_offset,'r--');
axis tight;
xlabel('log Amplitude');
legend('original','Calibrated');
subplot(2,2,4);
plot(phase,'k');
hold on
plot(femphase+phase_offset,'r--');
axis tight;
xlabel('Phase');
legend('original','Calibrated');

% restore NaNs to data
fem_datatmp = fem_data;
lnI_offsettmp = lnI_offset;
phase_offsettmp = phase_offset;
fem_data = [];
lnI_offset = [];
phase_offset = [];
datanum = 0;
for i = 1 : n
  for j = 1 : m
      if mesh.link(i,j) ~= 0
          datanum = datanum + 1;
          fem_data = [fem_data; fem_datatmp(datanum,:)];
          lnI_offset = [lnI_offset; lnI_offsettmp(datanum)];
          phase_offset = [phase_offset; phase_offsettmp(datanum)];
      else
          fem_data = [fem_data; NaN NaN];
          lnI_offset = [lnI_offset; NaN];
          phase_offset = [phase_offset; NaN];
      end
  end
end
