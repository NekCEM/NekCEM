
 clear all

 format long e
 load   'spatial_matrix_dtn1.dat';
 load   'eigenvalue_dtn1.dat';
 load   'spatial_matrix_dtn2.dat';
 load   'eigenvalue_dtn2.dat';

%-------------------------------------------
% Data copy and reshape
%-------------------------------------------

 xx20   = eigenvalue_dtn1(:,1);  % ying's grids for x
 yy20   = eigenvalue_dtn1(:,2);  % ying's grids for y
 xx21   = eigenvalue_dtn2(:,1);  % ying's grids for x
 yy21   = eigenvalue_dtn2(:,2);  % ying's grids for y

 u0    = spatial_matrix_dtn1(:);
 u1    = spatial_matrix_dtn2(:);

 
 
 neltx=3;
 nelty=2;
 nx1=4;
 ny1=4;
 
 %npts=length(x);
 npts=sqrt(length(u0))/2;
 
 npts2=npts*2;
 nn=npts2;

 uur0= reshape(u0,nn,nn);
 uur1= reshape(u1,nn,nn);
 uur20=uur0;
 uur21=uur1;
 
 pp=zeros(nelty*ny1,neltx*nx1);
 po=zeros(nelty*ny1,neltx*nx1);
 ind=0;
 
 for ky=1:nelty
     k2=(ky-1)*neltx*nx1*ny1;
     for kx=1:neltx
         k1=(kx-1)*nx1*ny1+k2;
         for jy=1:ny1
             j1=(jy-1)*nx1+k1
             for jx=1:nx1
               np=jx+j1
               po((ky-1)*ny1+jy,(kx-1)*nx1+jx)=np;
               
               if jx==nx1
                   
                   bb(ind+1)=np;
                   bb(ind+2)=np+npts;
                   ind=ind+2;
                   pp((ky-1)*ny1+jy,(kx-1)*nx1+jx)=1.0;
               elseif (jy==1 && ky>1)
                   
                   bb(ind+1)=np;
                   bb(ind+2)=np+npts;
                   ind=ind+2;
                   pp((ky-1)*ny1+jy,(kx-1)*nx1+jx)=1.0;
               end
            end
         end
     end    
 end
 bb=bb';
 bb=sort(bb);
 
 for k=ind:-1:1
    uur0(bb(k),:)=[];
    uur0(:,bb(k))=[];
    uur1(bb(k),:)=[];
    uur1(:,bb(k))=[];


 end


 cc=eig(uur1);
 xx1=real(cc);
 yy1=imag(cc);
 

 
 cond(uur1)
 

 
 
 
figure(300);set(gca,'fontsize',20);
box on
figure(300);hold on;spy(uur1,'b.');title('Spatial Operator: Structure');axis('square');
figure(300);print -depsc single_spy0103.eps;print -dpng single_spy0103.png

% draw eigen
ermax=max(xx1);eimax=max(yy1);
figure(400);set(gca,'fontsize',20);
figure(400);hold on;plot(xx1,yy1,'b.','MarkerSize', 16);title('Eigenvalue Distribution (\lambda)');
figure(400);hold on;

box on

axis('square');
axis([-35 5 -.5 .5]);
figure(400);xlabel(['max: real(\lambda)=',num2str(ermax),', imag(\lambda)=',num2str(eimax)]);
figure(400);print -depsc single_eig0103.eps;print -dpng single_eig0103.png







figure(100);set(gca,'fontsize',20);
box on
figure(100);hold on;spy(uur20,'b.');title('Spatial Operator: Structure');axis('square');
figure(100);print -depsc single_spy0003.eps;print -dpng single_spy0003.png


% draw eigen

ermax=max(xx20);eimax=max(yy20);
figure(200);set(gca,'fontsize',20);
figure(200);hold on;plot(xx20,yy20,'b.','MarkerSize', 16);title('Eigenvalue Distribution (\lambda)');
figure(200);hold on;

box on

axis('square');
axis([-35 5 -.5 .5]);
figure(200);xlabel(['max: real(\lambda)=',num2str(ermax),', imag(\lambda)=',num2str(eimax)]);
figure(200);print -depsc single_eig0003.eps;print -dpng single_eig0003.png


uu=uur21-uur20;
figure(500);set(gca,'fontsize',20);
box on
figure(500);hold on;spy(uu,'b.');title('Spatial Operator: Structure');axis('square');
figure(500);print -depsc single_spy_dtn.eps;print -dpng single_spy_dtn.png




