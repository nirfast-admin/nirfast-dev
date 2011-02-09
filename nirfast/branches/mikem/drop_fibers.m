function [link] = drop_fibers(link,fibers)
% [amesh hmesh] = drop_fibers(link,fibers)
% Drops fibers from input meshes in the vector fibers.
% Edit to drop wavelengths if needed.

for i = 1:length(fibers)
    j = fibers(i);
    a=link(:,1)==j;
    b=link(:,2)==j;
    link(a,3:end)=0;
    link(b,3:end)=0;
end

end

