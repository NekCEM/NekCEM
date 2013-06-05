
clear all

format long e
load 'fort.10';

x =fort(:,2);  % ying's grids for x
y =fort(:,3);  % ying's grids for y
u =fort(:,4);  % ying's field input
ur=fort(:,5);  % interpolation for real
ui=fort(:,6);  % interpolation for imag

nn = length(x);
n  = sqrt(nn);
xx = reshape(x ,n,n);
yy = reshape(y ,n,n);
uu = reshape(u ,n,n);
uur= reshape(ur,n,n);
uui= reshape(ui,n,n);

maxmin1=[max(u ), min(u)]
maxmin2=[max(ur), min(ur)]
maxmin3=[max(ui), min(ui)]

tmax1=maxmin1(1);
tmax2=maxmin1(2);
tmax3=maxmin2(1);
tmax4=maxmin2(2);
tmax5=maxmin3(1);
tmax6=maxmin3(2);

maxerr_real=max(u-ur)
maxerr_imag=max(u-ui)


%figure(1);
% subplot(3,1,1);plot3(x,y,u );caxis('auto'); xlabel('input'); 
% subplot(3,1,2);plot3(x,y,ur);caxis('auto'); xlabel('interp: real');
% subplot(3,1,3);plot3(x,y,ui);caxis('auto'); xlabel('interp: imag');
 figure(2); title('comparison: interpolations'); 
 subplot(3,1,1);mesh(uu ); xlabel(['input: max=',num2str(tmax1),' min=',num2str(tmax2)]); view(2) 
 subplot(3,1,2);mesh(uur); xlabel(['input: max=',num2str(tmax3),' min=',num2str(tmax4)]); view(2) 
 subplot(3,1,3);mesh(uui); xlabel(['input: max=',num2str(tmax5),' min=',num2str(tmax6)]); view(2) 


