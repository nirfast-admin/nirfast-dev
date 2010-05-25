function save_data(data,fn)

% save_data(data,fn)
%
% where data is the structured variable containing 
% amplitude and phase (*.paa)
% Saves calculated boundary data. Automatically determines type of data 




% check if the user reversed the inputs
if ischar(data)
    temp = data;
    data = fn;
    fn = temp;
end

% SPECTRAL
if isfield(data,'wv') == 1
  fid = fopen(fn,'w');
  [nrow,ncol]=size(data.wv);
  for i = 1 : nrow
    for j = 1 : ncol
      fprintf(fid,'w%-10g\t',data.wv(i,j));
    end
    fprintf(fid,'\n');
  end
  
  [nrow,ncol]=size(data.paa);
  for i = 1 : nrow
    for j = 1 : ncol
      fprintf(fid,'%-10g\t',data.paa(i,j));
    end
    fprintf(fid,'\n');
  end
  
  
% FLUORESCENCE
elseif isfield(data,'amplitudex') || isfield(data,'amplitudemm') || isfield(data,'amplitudefl')...
        || isfield(data,'phasex') || isfield(data,'phasemm') || isfield(data,'phasefl')
    fid = fopen(fn,'w');
    str = '';
    str2 = '';
    paa = [];
    
    % check which fields exist
    if isfield(data,'phasex')
        str = [str ',''xphase'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.phasex];
    end
    if isfield(data,'amplitudex')
        str = [str ',''xamplitude'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.amplitudex];
    end
    if isfield(data,'phasefl')
        str = [str ',''flphase'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.phasefl];
    end
    if isfield(data,'amplitudefl')
        str = [str ',''flamplitude'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.amplitudefl];
    end
    if isfield(data,'phasemm')
        str = [str ',''mmphase'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.phasemm];
    end
    if isfield(data,'amplitudemm')
        str = [str ',''mmamplitude'''];
        str2 = [str2 '%-10s\t'];
        paa = [paa data.amplitudemm];
    end
    
    % print file header
    str = ['fprintf(fid,''' str2 '\n''' str ');'];
    eval(str);
    
    % print data
    [nrow,ncol]=size(paa);
    for i = 1 : nrow
      for j = 1 : ncol
        fprintf(fid,'%-10g\t',paa(i,j));
      end
      fprintf(fid,'\n');
    end
  

% STANDARD
else
  fid = fopen(fn,'w');
 [nrow,ncol]=size(data.paa);
  for j = 1 : ceil(ncol/2)
    fprintf(fid,...
	    '%-10s\t%-10s\t',...
	    'amplitude','phase');
  end
  fprintf(fid,'\n');
  for i = 1 : nrow
    for j = 1 : ncol
      fprintf(fid,'%-10g\t',data.paa(i,j));
    end
    fprintf(fid,'\n');
  end
end
fclose(fid);
