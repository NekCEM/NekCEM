
load fort.60

nn =size(fort(:,1),1) 
n  =sqrt(nn(:,1)) 
a  =reshape(fort(:,1),n,n);
m  =n/3;

a11 = a(1+0*m:1*m,1+0*m:1*m);
a12 = a(1+0*m:1*m,1+1*m:2*m);
a13 = a(1+0*m:1*m,1+2*m:3*m);

a21 = a(1+1*m:2*m,1+0*m:1*m);
a22 = a(1+1*m:2*m,1+1*m:2*m);
a23 = a(1+1*m:2*m,1+2*m:3*m);

a31 = a(1+2*m:3*m,1+0*m:1*m);
a32 = a(1+2*m:3*m,1+1*m:2*m);
a33 = a(1+2*m:3*m,1+2*m:3*m);

figure(11); spy(a);

dd  = -1.5:0.01:1.5;
ee  =  eig(a);

figure(22)
set(gca,'fontsize',13);
plot(ee,'k.');axis('square');%grid; 
figure(22);hold on;plot(0,dd,'r--');axis('square');%grid; 
dd=max(real(ee));
xlabel(['max real eig=',num2str(dd)])


