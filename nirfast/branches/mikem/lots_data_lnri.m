
wv = amesh.wv;
for i = 1:2:12
plot_lnri(adata.paa(:,i:i+1),amesh,wv((i+1)/2));
end