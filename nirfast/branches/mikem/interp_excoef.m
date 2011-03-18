d = importdata('scottexcoef.txt');
clear exc
clf
plot(d(:,1),d(:,2),'bo')
hold on
plot(d(:,1),d(:,3),'ro')
plot(d(:,1),d(:,4),'go')

x = 650:950;
pp2 = spline(d(:,1),d(:,2));
pp3 = spline(d(:,1),d(:,3));
pp4 = spline(d(:,1),d(:,4));

% these start at 650
hbo = ppval(x,pp2);
hb = ppval(x,pp3);
water = ppval(x,pp4);

plot(x,hbo,'b')
plot(x,hb,'r')
plot(x,water,'g')

% wv = [661 735 785 808 826 849 903 912 948];
wv = [661 735 785 808 826 849];
exc(:,1) = hbo(wv - 649);
exc(:,2) = hb(wv - 649);
exc(:,3) = water(wv - 649);