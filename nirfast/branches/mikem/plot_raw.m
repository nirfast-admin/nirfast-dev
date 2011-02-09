function plot_raw(data,mesh,color)
% Plot Raw plots data as it comes out of the system with NaN's and 0's
% accounted for. Those are removed from the paa before plotting.

a = whos('data');
if strcmp(a.class,'struct')==1
    data = data.paa;
end

a = mesh.link(:)==0;
data(a,:) = NaN;
p = ~isnan(data(:,1));
data = data(p,:);


figure
clf
plot(log(data(:,1)),color);
title('amplitude')

figure
clf
plot(data(:,2),color);
title('phase')

end

