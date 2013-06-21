
 clear all

 format long e
 load   'fort.10';

%-------------------------------------------
% Data Sizes
%-------------------------------------------
 nn1 = fort(1,1); % size of grids in x
 nn2 = fort(1,2); % size of grids in y
 nn3 = fort(1,3); % size of grids in z
 nn  = fort(1,4); % length of field
 
%-------------------------------------------
% Data copy and reshape
%-------------------------------------------
 fst =2; lst=nn+1;
 x   = fort(fst:lst,2);  % ying's grids for x
 y   = fort(fst:lst,3);  % ying's grids for y
 z   = fort(fst:lst,4);  % ying's grids for z
 ur  = fort(fst:lst,5);  % interpolated field: real
 ui  = fort(fst:lst,6);  % interpolated field: imag
 ur0 = fort(fst:lst,7);  % ying's field input: real
 ui0 = fort(fst:lst,8);  % ying's field input: imag

 xx  = reshape(x  ,nn2,nn1);
 yy  = reshape(y  ,nn2,nn1);
 uur = reshape(ur ,nn2,nn1);
 uui = reshape(ui ,nn2,nn1);
 uur0= reshape(ur0,nn2,nn1);
 uui0= reshape(ui0,nn2,nn1);
 
%-------------------------------------------
% Get max/min values of the fields
%-------------------------------------------
 maxmin_ur0=[max(ur0), min(ur0)];
 maxmin_ui0=[max(ui0), min(ui0)];
 maxmin_ur =[max(ur ), min(ur )];
 maxmin_ui =[max(ui ), min(ui )];

 ur0max = maxmin_ur0(1);
 ur0min = maxmin_ur0(2);
 ui0max = maxmin_ui0(1);
 ui0min = maxmin_ui0(2);

 urmax  = maxmin_ur (1);
 urmin  = maxmin_ur (2);
 uimax  = maxmin_ui (1);
 uimin  = maxmin_ui (2);

 max_min_ur0=[ur0max ur0min] 
 max_min_ur =[urmax urmin] 

 max_min_ui0=[ui0max ui0min] 
 max_min_ui =[uimax uimin] 

% Compute pointwise errors of the fields
%-------------------------------------------

 err_ur  = abs(uur0-uur);
 err_ui  = abs(uui0-uui);

 maxerr_ur= norm(err_ur,Inf) ;
 maxerr_ui= norm(err_ui,Inf) ;
 l2err_ur = norm(err_ur,  2) ;
 l2err_ui = norm(err_ui,  2) ;

 maxerr_real_imag=[maxerr_ur maxerr_ui]
 L2err_real_imag =[l2err_ur  l2err_ui ]

%-------------------------------------------
% Draw Figures
%-------------------------------------------
 xmax= max(max(xx)); xmin= min(min(xx));  % get range in x
 ymax= max(max(yy)); ymin= min(min(yy));  % get range in y

%-------------------------------------------
% Use contour:: Imaginary Part: fields and pointwise errors
%-------------------------------------------
 figure(1)     ; title('Imag Part'); 
 subplot(3,1,1); contour(xx,yy,uui0,50); xlabel(['Imag:: TFE: max=',num2str(ui0max),'; min=',num2str(ui0min)]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); contour(xx,yy,uui,50); xlabel(['Imag:: SEM: max=',num2str(uimax),'; min=',num2str(uimin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); contour(xx,yy,err_ui,50); xlabel(['Imag:: Pointwise Errors: max|TFE-SEM|=',num2str(maxerr_ui)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 figure(1); print -depsc imag.eps
 figure(1); print -dpng  imag.png

%-------------------------------------------
% Use contour:: Real Part: fields and pointwise errors
%-------------------------------------------
 figure(2)     ; title('Real Part'); 
 subplot(3,1,1); contour(xx,yy,uur0,50); xlabel(['Real:: TFE: max=',num2str(ur0max),'; min=',num2str(ur0min)]); 
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,2); contour(xx,yy,uur,50); xlabel(['Real:: SEM: max=',num2str(urmax),'; min=',num2str(urmin)]); 
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 subplot(3,1,3); contour(xx,yy,err_ur,50); xlabel(['Real:: Pointwise Errors: max|TFE-SEM|=',num2str(maxerr_ur)]); 
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y') 

 figure(2); print -depsc real.eps
 figure(2); print -dpng  real.png

%-------------------------------------------
% Use mesh:: Imaginary Part: fields and pointwise errors
%-------------------------------------------
 figure(3)     ; title('Imag Part');
 subplot(3,1,1); mesh(xx,yy,uui0); xlabel(['Imag:: TFE: max=',num2str(ui0max),'; min=',num2str(ui0min)]);
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 subplot(3,1,2); mesh(xx,yy,uui); xlabel(['Imag:: SEM: max=',num2str(uimax),'; min=',num2str(uimin)]);
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 subplot(3,1,3); mesh(xx,yy,err_ui); xlabel(['Imag:: Pointwise Errors: max|TFE-SEM|=',num2str(maxerr_ui)]);
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 figure(3); print -depsc imag2.eps
 figure(3); print -dpng  imag2.png

%-------------------------------------------
% Use mesh:: Real Part: fields and pointwise errors
%-------------------------------------------
 figure(4)     ; title('Real Part');
 subplot(3,1,1); mesh(xx,yy,uur0); xlabel(['Real:: TFE: max=',num2str(ur0max),'; min=',num2str(ur0min)]);
 subplot(3,1,1); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 subplot(3,1,2); mesh(xx,yy,uur); xlabel(['Real:: SEM: max=',num2str(urmax),'; min=',num2str(urmin)]);
 subplot(3,1,2); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 subplot(3,1,3); mesh(xx,yy,err_ur); xlabel(['Real:: Pointwise Errors: max|TFE-SEM|=',num2str(maxerr_ur)]);
 subplot(3,1,3); view(2); axis([xmin xmax ymin ymax]); colorbar; ylabel('y')

 figure(4); print -depsc real2.eps
 figure(4); print -dpng  real2.png
%-------------------------------------------
