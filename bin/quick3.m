
load fort.10

nn =size(fort(:,1),1) 
n  =sqrt(nn(:,1)) 
a  =reshape(fort(:,1),n,n);

ee  =  eig(a);

figure(22)
set(gca,'fontsize',13);
plot(ee,'k.');axis('square');%grid; 


