
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

tmax1=maxmin1(1);
tmax2=maxmin1(2);
tmax3=maxmin2(1);
tmax4=maxmin2(2);
tmax5=maxmin3(1);
tmax6=maxmin3(2);

err_real= (uu-uur);
err_imag= (uu-uui);

errmax1= max(max(uu-uur))
errmin1= min(min(uu-uur))
errmax2= max(max(uu-uui))
errmin2= min(min(uu-uui))

xmax=max(max(xx)); xmin=min(min(xx));
ymax=max(max(yy)); ymin=min(min(yy));

%figure(1);
% subplot(3,1,1);plot3(x,y,u );caxis('auto'); xlabel('input'); 
% subplot(3,1,2);plot3(x,y,ur);caxis('auto'); xlabel('interp: real');
% subplot(3,1,3);plot3(x,y,ui);caxis('auto'); xlabel('interp: imag');

 figure(2); title('comparison: interpolations'); 
 subplot(3,1,1);mesh(xx,yy,uu ); xlabel(['imag   input: max=',num2str(tmax1),'; min=',num2str(tmax2)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 
 subplot(3,1,2);mesh(xx,yy,uui); xlabel(['imag numeric: max=',num2str(tmax5),'; min=',num2str(tmax6)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 
 subplot(3,1,3);mesh(xx,yy,uur); xlabel(['real numeric: max=',num2str(tmax3),'; min=',num2str(tmax4)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 figure(3); title('Error Distributions: all are imaginary part');
 subplot(3,1,1);mesh(xx,yy,err_imag); xlabel(['err imag: max=',num2str(errmax2),'; min=',num2str(errmin2)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 
 subplot(3,1,2);mesh(xx,yy,err_imag); xlabel(['err imag: max=',num2str(errmax2),'; min=',num2str(errmin2)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 
 subplot(3,1,3);mesh(xx,yy,err_imag); xlabel(['err imag: max=',num2str(errmax2),'; min=',num2str(errmin2)]); 
 view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 


