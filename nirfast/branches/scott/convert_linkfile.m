function convert_linkfile(fn_link, fn_newlink)

link = load(fn_link);
[n,m] = size(link);

if nargin == 2
    fn_link = fn_newlink;
end

fid = fopen(fn_link,'w');
fprintf(fid,'%s\n','source det active');
fl = 0;
for i = 1:n
    for j = 1:m
        if link(i,j) ~= 0
            fprintf(fid, '%g %g %g\n', i, link(i,j), 1);
        elseif link(i,j) == 0
            fprintf(fid, '%g %g %g\n', i, NaN, 0);
            fl = 1;
        end
    end
end
fclose(fid);

if fl == 1
disp('WARNING: One or more element values = zero in the original link file!')
disp('convert_linkfile.m could not determine what detector # was set to zero.')
disp('The new link file was written with an NaN in element.')
disp('Please view and adjust manually.')
end
