function [reg proj_error] = read_logfile(fn)
fid = fopen(fn,'r');
it=1;
while 1
    tline = fgetl(fid);
    log{it,1} = tline;
    if ~ischar(tline),   break,   end
    it = it+1;
end
fclose(fid);

sol = str2num(log{it-1});
for i = 1:it
    if strcmp(['Iteration Number          = ',num2str(sol)],log{i})==1;
       line = log{i+1}; 
       break
    end
end

proj_error = line(29:end);
reg = str2num(log{3}(22:end));
