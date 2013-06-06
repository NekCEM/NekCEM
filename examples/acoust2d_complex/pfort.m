
 clear all

 format long e
 load 'fort.10';

%-------------------------------------------
% Data Sizes
%-------------------------------------------
 nn1= fort(1,1); % size of grids in x
 nn2= fort(1,2); % size of grids in y
 nn3= fort(1,3); % size of grids in z
 nn = fort(1,4); % length of field
 
%-------------------------------------------
% Data copy and reshape
%-------------------------------------------
 x  = fort(2:nn+1,2);  % ying's grids for x
 y  = fort(2:nn+1,3);  % ying's grids for y
 z  = fort(2:nn+1,4);  % ying's grids for z
 ur = fort(2:nn+1,5);  % interpolated field: real
 ui = fort(2:nn+1,6);  % interpolated field: imag
 ur0= fort(2:nn+1,7);  % ying's field input: real
 ui0= fort(2:nn+1,8);  % ying's field input: imag

 xx  = reshape(x  ,nn2,nn1);
 yy  = reshape(y  ,nn2,nn1);
 uur = reshape(ur ,nn2,nn1);
 uui = reshape(ui ,nn2,nn1);
 uur0= reshape(ur0,nn2,nn1);
 uui0= reshape(ui0,nn2,nn1);

%-------------------------------------------
% Get max/min values of the fields
%-------------------------------------------
 maxmin_ur0=[max(ur0), min(ur0)]
 maxmin_ui0=[max(ui0), min(ui0)]
 maxmin_ur =[max(ur ), min(ur )]
 maxmin_ui =[max(ui ), min(ui )]

 ur0max = maxmin_ur0(1);
 ur0min = maxmin_ur0(2);
 ui0max = maxmin_ui0(1);
 ui0min = maxmin_ui0(2);

 urmax  = maxmin_ur (1);
 urmin  = maxmin_ur (2);
 uimax  = maxmin_ui (1);
 uimin  = maxmin_ui (2);

%-------------------------------------------
% Compute pointwise errors of the fields
%-------------------------------------------
 err_ur = abs((uur0-uur));
 err_ui = abs((uui0-uui));

 errmax_ur= abs(max(max(uur0-uur)))
 errmin_ur= abs(min(min(uur0-uur)))
 errmax_ui= abs(max(max(uui0-uui)))
 errmin_ui= abs(min(min(uui0-uui)))

%-------------------------------------------
% Draw Figures
%-------------------------------------------
 xmax= max(max(xx)); xmin= min(min(xx));  % get range in x
 ymax= max(max(yy)); ymin= min(min(yy));  % get range in y

%-------------------------------------------
% Imaginary Part: fields and pointwise errors
%-------------------------------------------
 figure(1); title('Imag Part'); 
 subplot(3,1,1); contour(xx,yy,uui0,40); xlabel(['Imag:: FTE: max=',num2str(ui0max),'; min=',num2str(ui0min)]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); contour(xx,yy,uui,40); xlabel(['Imag:: SEM: max=',num2str(uimax),'; min=',num2str(uimin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); contour(xx,yy,err_ui,40); xlabel(['Imag:: Pointwise Errors: max|FTE-SEM|=',num2str(errmax_ui)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 figure(1); print -depsc imag.eps
 figure(1); print -dpng  imag.png
%-------------------------------------------
% Real Part: fields and pointwise errors
%-------------------------------------------
 figure(2); title('Real Part'); 
 subplot(3,1,1); contour(xx,yy,uur0,40); xlabel(['Real:: FTE: max=',num2str(ur0max),'; min=',num2str(ur0min)]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); contour(xx,yy,uur,40); xlabel(['Real:: SEM: max=',num2str(urmax),'; min=',num2str(urmin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); contour(xx,yy,err_ur,40); xlabel(['Real:: Pointwise Errors: max|FTE-SEM|=',num2str(errmax_ur)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 figure(2); print -depsc real.eps
 figure(2); print -dpng  real.png
%-------------------------------------------
