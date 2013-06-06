
 clear all

 format long e
 load 'fort.10';

%-------------------------------------------
% Data Sizes
%-------------------------------------------
 nn1= fort(1,1); % size of grids in x
 nn2= fort(1,2); % size of grids in y
 nn = fort(1,3);
 
%-------------------------------------------
% Data copy and reshape
%-------------------------------------------
 x  = fort(2:nn+1,2);  % ying's grids for x
 y  = fort(2:nn+1,3);  % ying's grids for y
 u  = fort(2:nn+1,4);  % ying's field input
 ur = fort(2:nn+1,5);  % interpolation for real
 ui = fort(2:nn+1,6);  % interpolation for imag

 xx = reshape(x ,nn2,nn1);
 yy = reshape(y ,nn2,nn1);
 uu = reshape(u ,nn2,nn1);
 uur= reshape(ur,nn2,nn1);
 uui= reshape(ui,nn2,nn1);

%-------------------------------------------
% Get max/min values of the fields
%-------------------------------------------
 maxmin1=[max(u ), min(u )]
 maxmin2=[max(ur), min(ur)]
 maxmin3=[max(ui), min(ui)]

 umax   = maxmin1(1);
 umin   = maxmin1(2);
 urmax  = maxmin2(1);
 urmin  = maxmin2(2);
 uimax  = maxmin3(1);
 uimin  = maxmin3(2);

%-------------------------------------------
% Compute pointwise errors of the fields
%-------------------------------------------
 err_ur = (uu-uur)  ;
 err_ui = (uu-uui)  ;

 errmax_ur= max(max(uu-uur))
 errmin_ur= min(min(uu-uur))
 errmax_ui= max(max(uu-uui))
 errmin_ui= min(min(uu-uui))

%-------------------------------------------
% Draw Figures
%-------------------------------------------
 xmax= max(max(xx)); xmin= min(min(xx));  % get range in x
 ymax= max(max(yy)); ymin= min(min(yy));  % get range in y

%-------------------------------------------
% Imaginary Part: fields and pointwise errors
%-------------------------------------------
 figure(1); title('Imag Part'); 
 subplot(3,1,1); mesh(xx,yy,uu ); xlabel(['FTE: max=',num2str(umax ),'; min=',num2str(umin )]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); mesh(xx,yy,uui); xlabel(['SEM: max=',num2str(uimax),'; min=',num2str(uimin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); mesh(xx,yy,err_ui); xlabel(['Pointwise Errors: max|FTE-SEM|=',num2str(errmax_ui)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

%-------------------------------------------
% Real Part: fields and pointwise errors
%-------------------------------------------
 figure(2); title('Real Part'); 
 subplot(3,1,1); mesh(xx,yy,uu ); xlabel(['FTE: max=',num2str(umax ),'; min=',num2str(umin )]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); mesh(xx,yy,uur); xlabel(['SEM: max=',num2str(urmax),'; min=',num2str(urmin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); mesh(xx,yy,err_ur); xlabel(['Pointwise Errors: max|FTE-SEM|=',num2str(errmax_ur)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 
%-------------------------------------------
