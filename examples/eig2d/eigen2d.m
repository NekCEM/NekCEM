load eig2d_cen.dat
load eig2d_up.dat
load str2d_cen.dat
load str2d_up.dat

% construct spatial matrix "m":
s   = str2d_cen(:);
slen= length(s);slen=sqrt(slen);
m   =  reshape(s,slen,slen);
figure(100);set(gca,'fontsize',18);
figure(100);hold on;spy(m,'k.');title('2D: Structure of Spatial Operator (Central)');axis('square');
figure(100);print -depsc str2d_cen.eps;print -dpng str2d_cen.png

% construct spatial matrix "m":
s   = str2d_up(:);
slen= length(s);slen=sqrt(slen);
m   =  reshape(s,slen,slen);
figure(200);set(gca,'fontsize',18);
figure(200);hold on;spy(m,'k.');title('2D: Structure of Spatial Operator (Upwind)');axis('square')
figure(200);print -depsc str2d_up.eps;print -dpng str2d_up.png

% draw eigen
e=eig2d_cen(:,:); ermax=max(e(:,1));eimax=max(e(:,2));
figure(300);set(gca,'fontsize',18);
figure(300);hold on;plot(e(:,1),e(:,2),'k.');title('2D: Eigenvalue Distribution (Central)');
figure(300);hold on;axis('square');axis([-1 1 -1.0 1.0]);
figure(300);xlabel(['max: real(\lambda)=',num2str(ermax),', imag(\lambda)=',num2str(eimax)]);
figure(300);print -depsc eig2d_cen.eps;print -dpng eig2d_cen.png

e=eig2d_up(:,:); ermax=max(e(:,1));eimax=max(e(:,2));
figure(400);set(gca,'fontsize',18);
figure(400);hold on;plot(e(:,1),e(:,2),'k.');title('2D: Eigenvalue Distribution (Upwind)');
figure(400);hold on;axis('square');axis([-1 1 -1.0 1.0]);
figure(400);xlabel(['max: real(\lambda)=',num2str(ermax),', imag(\lambda)=',num2str(eimax)]);
figure(400);print -depsc eig2d_up.eps;print -dpng eig2d_up.png



