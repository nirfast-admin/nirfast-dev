function convert2spectraldata_mike(fn,op_fn,rep,wv_array)
% convert2spectraldata
% Converts paa files from nir/mr system to spectral paa files.
% Give (fn,op_fn,rep,wv_array)
% fn is the base name of the calibrated data (usually '..._op.paa')
%  if fn is given as a matrix of data, [240 x #wv*2], it will just add header and save. 
% op_fn will be the converted data basename (e.g. 'op_fn' -> 'op_fn.paa')
% Rep is the repitition number
% wv_array is the wvlengths to use (enter as doubles)- default is [661 761 785 808 826 849]

if(nargin <4)
    wv_array = [661 761 785 808 826 849];
end

% Load data if it's a string
if ischar(fn) == 1
    data_op = [];
    for(i=1:length(wv_array))
        data = load([fn num2str(wv_array(i)) 'nm_rep' num2str(rep) '.paa']);
        data_op = [data_op data(:,1) data(:,2)];
    end
else 
    data_op = fn;
end


% Make header
hdr_op = [];
for(i=1:length(wv_array))
    hdr_op = [hdr_op 'w' num2str(wv_array(i)) ' '];
end


% Write new data files    
[nrow,ncol]=size(data_op);
fid = fopen([op_fn 'rep' num2str(rep) '.paa'],'w');
fprintf(fid,'%s\n',hdr_op);

for i = 1 : nrow
  for j = 1 : ncol
    fprintf(fid,'%g ',data_op(i,j));
  end
      fprintf(fid,'\n');
end
fclose(fid);