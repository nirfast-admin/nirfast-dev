% converts new meshes and data to old format

function [data mesh] = new2old_link(data, mesh)

[m,n] = size(data.link); 
j=1;

for i = 3:n
    a = data.link(:,i) == 0;
    data.paa(a,j:j+1) = NaN;
    j = j+2;
end

for i = 1:16
a = mesh.link(:,1) == i;
link(i,:) = mesh.link(a,2);
end

mesh.link = link;
