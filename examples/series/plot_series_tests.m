data=load('series_mean.dat');
[m,n]=size(data);

figure(2); clf;
plot(data(:,1),data(:,2),'k');

data=load('series_distr.dat');
[m,n]=size(data);

figure(1); clf;
plot(data(:,1),data(:,2),'k');

data_diff=load('series_diff.dat');

cmpmat=load('series_cmp.dat');

figure(4); clf;
contourf(cmpmat);
set(gca,'FontSize',24);
colormap(hot);
cmap=flipud(colormap);
colormap(cmap);