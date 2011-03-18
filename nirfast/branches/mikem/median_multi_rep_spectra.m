function median_multi_rep_spectra(data_fn,save_fn,num_reps,source)

% Load a series of repetitions and determine the median at each pixel.
% Save result as save_fn


for i = 1:num_reps
    % load basis data for one rep
    data_temp = load([data_fn,sprintf('_s%g_rep%g.raw',source,i)]);
    [n,m] = size(data_temp);
    data(:,:,i) = data_temp;
end

[c,r,d] = size(data);
data_avg = median(data,3);
data_dev = std(data,1,3);

save([save_fn,'.raw'],'data_avg','-ASCII');
copyfile([data_fn,'.log'],[save_fn,'.log']);