function param = input_calib

param.fl_baseline_region_pixels = [20:220];

% 
% 
% param.fl_searchpeak_range_AF750 = [520:600];
% param.fl_searchpeak_range_Licor800 = [620:660];
% 
% param.fl_peak_width_pixels = 100;
% param.fl_baseline_region_pixels = [200:220];
% 
% param.x_searchpeak_range = [80:200];
% 
% param.x_peak_width_pixels = 80;
% 
% param.Basis_autofluor = 'Basis_autofluor.txt';
% param.Basis_AF750 = 'Basis_AF750.txt';
% param.Basis_Licor800 = 'Basis_Licor.txt'
% param.specfitmin = [350];
%  
% % one spectrum for each spectrograph which provides the ratio of
% % intensities with and without the 720 nm filter.  Note data before pixel
% % #360 has been set = 1 (outside of filtering range), so do not use this
% % range in spectral unmixing
% param.flfilter_scalefactor = 'Filter720_scalefactor.txt';
% 
% % 16 values: inter-detector factors on fl/trans for AF750
% param.flx_ratio_scalefactor_drug1 = 'AF750_ratio_scalefactor.txt'; 
% % 16 values: inter-detector factors on fl/trans for Licor800
% param.flx_ratio_scalefactor_drug2 = 'Licor800_ratio_scalefactor.txt'; 
% 
% % Calculated value of the ND filter on the source, used for acquiring
% % excitation measurements.
% param.x_filter_OD = 4.7772;
% 
% 
% % FOR FUTURE USE?:
% % 16 values:  Pixel of fluor peak for each spectrograph 
% param.fl_peak_pixel_drug1 = '.txt';
% % 16 values:  Pixel of fluor peak for each spectrograph
% param.fl_peak_pixel_drug2 = '.txt';  
% % 16 values:  Pixel of laser peak for each spec
% param.x_peaks_pixel = '.txt';  