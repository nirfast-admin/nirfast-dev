
function [mua_init,mus_init,bem_data] = fit_data_slab(fn,...
						  data,...
						  frequency,...
						  iteration,link)

% fits data to a given model to give initial guesses for
% reconstruction as well as data offsets
% h dehghani 16 Jan 2002 hamid.dehghani@dartmouth.edu
%%modified for bem, subha 3/27/07
%%modified for slab geometry by Subha 11/21/07
[elem,nodes,source,param,meas,junk] = get_input_data(fn);

param_mus = (1./(3*param(:,3))) - param(:,2);

% convert log amplitude into amplitude
data(:,1) = exp(data(:,1));

% calculate the source / detector distance for each measurement
k = 1;r = 1; t = 1;
for i = 1 : length(source)
 for j = 1 : length(link(i,:))
   if link(i,j) ~= 0
     jj = link(i,j);
     %%separating reflectance measurements
     if (source(i,3) == meas(jj,3))
         dist_r(r,1) = sqrt(sum((source(i,:) - meas(jj,:)).^2));
         data_r(r,:) = data(k,:);
         index_meas(k) = 1;
         r = r+1;
     else
         %%separating transmittance measurements
         dist_t(t,1) = sqrt(sum((source(i,:) - meas(jj,:)).^2));
         data_t(t,:) = data(k,:);
         index_meas(k) = 0;
         t = t+1;
     end
     dist(k,1) = sqrt(sum((source(i,:) - meas(jj,:)).^2));
       k = k+1;
   end
 end
end

% Set lnrI, lnr and phase!
lnrI_ref = log(data_r(:,1).*dist_r);
lnrI_trans = log(data_t(:,1).*dist_t);
lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));
phase = data(:,2);

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
subplot(1,2,1);
plot(dist,lnrI,'.')
ylabel('lnrI');
xlabel('Source / Detector distance');
subplot(1,2,2);
plot(dist,phase,'.')
ylabel('Phase');
xlabel('Source / Detector distance');
drawnow
pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
m0 = polyfit(dist,phase,1); m0 = m0(1);
%%calculating slopes separately for reflectance and transmittance
%%measurements
m1_r = polyfit(dist_r,lnrI_ref,1); m1_r = m1_r(1);
m1_t = polyfit(dist_t,lnrI_trans,1); m1_t = m1_t(1);
%m1 = [m1_r m1_t]
%m1 = mean(m1);
m1 = m1_r;
disp('slope from reflectance = ');
m1_r
%disp('slope from transmittance = ');
%m1_t
m1_old = polyfit(dist,lnrI,1);m1_old = m1_old(1)
% fit data using an analytical model
% based on Pogue paper
omega = 2*pi*frequency*1e6;
ref_index = 1.33
c=(3e11./ref_index);

mua = 0.01; mus = 1; kappa = 1/(3*(mua+mus));

for i = 1 : 25
 f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
 alpha = -f1*cos(0.5*atan2(omega/c,mua));
 phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
 %%perturb mua and calculate slopes
 mua = mua + 0.0001; 
 kappa = 1/(3*(mua+mus));
 f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
 alpha1 = -f1*cos(0.5*atan2(omega/c,mua));
 phi1 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
 %%unperturb mua
 mua = mua - 0.0001;   
 %%calculate update to mua
 da = (alpha-alpha1)/0.0001;
 ds = (phi-phi1)/0.0001;
 mua = mua - (m1-alpha)/da;
 mus = mus - (m0-phi)/ds*0.1;
 kappa = 1/(3*(mua+mus));
 %%calculate the new slopes
 f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
 alpha = -f1*cos(0.5*atan2(omega/c,mua));
 phi = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
 %%perturb mus and calculate slopes
 mus = mus + 0.001; 
 kappa = 1/(3*(mua+mus));
 f1 = ((mua^2 + (omega/c)^2) / (kappa^2))^0.25;
 alpha2 = -f1*cos(0.5*atan2(omega/c,mua));
 phi2 = f1*sin(0.5*atan2(omega/c,mua))*180/pi;
 %%unperturb mus (back to original value)
 mus = mus - 0.001; 
 %%calculate updates
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
param(:,2) = mua;
param_mus(:) = mus;
param(:,3) = kappa;
%calib_mesh = strcat(fn,'_calib_mesh');
%save_mesh_bem(elem,nodes,param,source,meas,link,calib_mesh);
mua_init = mua;
mus_init = mus;
kappa_init = kappa;
%%save_mesh_bem(elem,nodes,param,source,meas,link,recon_mesh);


% Fit for mua and mus using FEM
jj = 0;
 ind_r = find(index_meas(:) == 1);
 ind_t = find(index_meas(:) == 0);
 disp('numerical fit ................');

while jj ~= iteration
 [junk,bem_data]=bem_dot_3d(fn,frequency,param);
 %fem_data = fem_data.paa;

 bemlnrI_r = log(bem_data(ind_r,1).*dist_r);
 bemlnrI_t = log(bem_data(ind_t,1).*dist_t);  
 bemphase = bem_data(:,2);
 %%calculate slopes
 phi0 = polyfit(dist,bemphase,1); phi0 = phi0(1);
 alpha0_r = polyfit(dist_r,bemlnrI_r,1); alpha0_r = alpha0_r(1);
 alpha0_t = polyfit(dist_t,bemlnrI_t,1); alpha0_t = alpha0_t(1);
 %alpha0 = [alpha0_r alpha0_t];
 %alpha0 = mean(alpha0);
 alpha0 = alpha0_r;
 disp('slope from reflectance = ');
   alpha0_r
 %disp('slope from transmittance = ');
 %  alpha0_t
 %%perturb mua and calculate slopes
param(:,2) = param(:,2) + 0.0001;
param(:,3) = 1./(3*(param(:,2)+param_mus));

 [junk,bem_data]=bem_dot_3d(fn,frequency,param);

 bemlnrI_r = log(bem_data(ind_r,1).*dist_r);
 bemlnrI_t = log(bem_data(ind_t,1).*dist_t);  
 bemphase = bem_data(:,2);

 phi1 = polyfit(dist,bemphase,1); phi1 = phi1(1);
 alpha1_r = polyfit(dist_r,bemlnrI_r,1); alpha1_r = alpha1_r(1);
 alpha1_t = polyfit(dist_t,bemlnrI_t,1); alpha1_t = alpha1_t(1);
 %alpha1 = [alpha1_r alpha1_t];
 %alpha1 = mean(alpha1);
 alpha1 = alpha1_r;
 %alpha1 = polyfit(dist,bemlnrI,1); alpha1 = alpha1(1);
 %%unperturb mua back to original value
 param(:,2) = param(:,2)-0.0001;

 %%calculate updates
 da = (alpha0-alpha1)/0.0001;
 ds = (phi0-phi1)/0.0001;

 param(:,2) = param(:,2) - (m1-alpha0)/da;
 param_mus = param_mus - (m0-phi0)/ds.*0.1;
param(:,3) = 1./(3*(param(:,2)+param_mus));
[param(1,2) param_mus(1)];
%%calculate the new slopes 
[junk,bem_data]=bem_dot_3d(fn,frequency,param);

 bemlnrI_r = log(bem_data(ind_r,1).*dist_r);
 bemlnrI_t = log(bem_data(ind_t,1).*dist_t);  
 bemphase = bem_data(:,2);

 phi0 = polyfit(dist,bemphase,1); phi0 = phi0(1);
 alpha0_r = polyfit(dist_r,bemlnrI_r,1); alpha0_r = alpha0_r(1);
 alpha0_t = polyfit(dist_t,bemlnrI_t,1); alpha0_t = alpha0_t(1);
 %alpha0 = [alpha0_r alpha0_t];
 %alpha0 = mean(alpha0);
 alpha0 = alpha0_r;
 %%perturb mus and calculate slopes
 param_mus(:) = param_mus(:)+0.001;
 param(:,3) = 1./(3*(param(:,2)+param_mus));

[junk,bem_data]=bem_dot_3d(fn,frequency,param);

 bemlnrI_r = log(bem_data(ind_r,1).*dist_r);
 bemlnrI_t = log(bem_data(ind_t,1).*dist_t);  
 bemphase = bem_data(:,2);

 phi2 = polyfit(dist,bemphase,1); phi2 = phi2(1);
 alpha2_r = polyfit(dist_r,bemlnrI_r,1); alpha2_r = alpha2_r(1);
 alpha2_t = polyfit(dist_t,bemlnrI_t,1); alpha2_t = alpha2_t(1);
 %alpha2 = [alpha2_r alpha2_t];
 %alpha2 = mean(alpha2);
 alpha2 = alpha2_r;
 %alpha2 = polyfit(dist,bemlnrI,1); alpha2 = alpha2(1);

 %%unperturb mus back to original value
 param_mus(:) = param_mus(:)-0.001;

%%calculate updates  
 da = (alpha0-alpha2)/0.001;
 ds = (phi0-phi2)/0.001;

 param(:,2) = param(:,2) - (m1-alpha0)/da*0.01;
 %param(find(param(:,2)<0),2) = param(find(param(:,2)<0),2) + da;
 param_mus(:) = param_mus(:) - (m0-phi0)/ds.*0.2;
  param(:,3) = 1./(3*(param(:,2)+param_mus));

 err_a = abs(mean(param(:,2))-mua_init)./mua_init;
 err_s = abs(mean(param_mus)-mus_init)./mus_init;

 if ((err_a < 0.001) & (err_s < 0.001))
   jj = iteration;
 else
   jj = jj + 1;
 end

 mua_init = mean(param(:,2));
 mus_init = mean(param_mus);
 disp('Global values calculated from Numerical fit');
 disp(['Iteration ' num2str(jj) ' of ' num2str(iteration)]);
 disp(['Absorption = ' num2str(mua_init) ' mm-1 with error of ' num2str(err_a)]);
 disp(['Scatter    = ' num2str(mus_init) ' mm-1 with error of ' num2str(err_s)]);
 disp('-------------------------------------------------');
end
