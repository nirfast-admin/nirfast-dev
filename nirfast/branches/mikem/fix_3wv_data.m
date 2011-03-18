% find extra points in link
a = homog_data_cw.link(:,1)~=homog_data_cw.link(:,2);
% homog
data.link = homog_data_cw.link(a,:);
data.paa = homog_data_cw.paa(a,1:2:5);
homog_data_cw = data;
homog_data_cw.wv = [903 912 948];
% anom
data.paa = anom_data_cw.paa(a,1:2:5);
data.link = anom_data_cw.link(a,:);
anom_data_cw = data;
anom_data_cw.wv = [903 912 948];
% combine homog
clear data
data.wv = [hdata2.wv homog_data_cw.wv];
data.link = [hdata2.link homog_data_cw.link(:,3:end)];
data.paa = [hdata2.paa homog_data_cw.paa];
homog_data_big = data;
% combine anom
clear data
data.wv = [adata2.wv anom_data_cw.wv];
data.link = [adata2.link anom_data_cw.link(:,3:end)];
data.paa = [adata2.paa anom_data_cw.paa];
anom_data_big = data;