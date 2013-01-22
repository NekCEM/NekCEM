
load spatial_matrix.dat
load eigenvalue.dat


% construct spatial matrix "m":
s   = spatial_matrix(:);
slen= length(s);slen=sqrt(slen);
m=  reshape(s,slen,slen); 
figure(100);spy(m);title('Spatial Operator: Strucutre');

e=eigenvalue(:,:); ermax=max(e(:,1));eimax=max(e(:,2));
figure(200);plot(e(:,1),e(:,2),'.');title('Eigenvalue Distribution (\lambda)');    
figure(200);xlabel(['max: real(\lambda)=',num2str(ermax),', imag(\lambda)=',num2str(eimax)]);
