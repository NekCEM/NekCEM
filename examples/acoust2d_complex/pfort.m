
clear all

format long e
load 'fort.10';

nn1=fort(1,1); % size of grids in x
nn2=fort(1,2); % size of grids in y
nn =fort(1,3);

x  =fort(2:nn+1,2);  % ying's grids for x
y  =fort(2:nn+1,3);  % ying's grids for y
u  =fort(2:nn+1,4);  % ying's field input
ur =fort(2:nn+1,5);  % interpolation for real
ui =fort(2:nn+1,6);  % interpolation for imag

xx = reshape(x ,nn2,nn1);
yy = reshape(y ,nn2,nn1);
uu = reshape(u ,nn2,nn1);
uur= reshape(ur,nn2,nn1);
uui= reshape(ui,nn2,nn1);

maxmin1=[max(u ), min(u)]
maxmin2=[max(ur), min(ur)]
maxmin3=[max(ui), min(ui)]

maxerr_real=max(u-ur)
maxerr_imag=max(u-ui)

tmax1=maxmin1(1);
tmax2=maxmin1(2);
tmax3=maxmin2(1);
tmax4=maxmin2(2);
tmax5=maxmin3(1);
tmax6=maxmin3(2);


%figure(1);
% subplot(3,1,1);plot3(x,y,u );caxis('auto'); xlabel('input'); 
% subplot(3,1,2);plot3(x,y,ur);caxis('auto'); xlabel('interp: real');
% subplot(3,1,3);plot3(x,y,ui);caxis('auto'); xlabel('interp: imag');
 figure(2); title('comparison: interpolations'); 
 subplot(3,1,1);mesh(uu ); xlabel(['input: max=',num2str(tmax1),' min=',num2str(tmax2)]); view(2) 
 subplot(3,1,2);mesh(uur); xlabel(['input: max=',num2str(tmax3),' min=',num2str(tmax4)]); view(2) 
 subplot(3,1,3);mesh(uui); xlabel(['input: max=',num2str(tmax5),' min=',num2str(tmax6)]); view(2) 


