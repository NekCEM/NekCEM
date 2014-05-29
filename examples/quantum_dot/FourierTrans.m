
load fourier.dat
load pulse.dat   

f1=fourier;
t2=pulse  ;
w = f1(:,1)
Fr= f1(:,2)
Fi= f1(:,3)
t = t2(:,1)
p = t2(:,2)

figure(100);plot(w,Fr,w,Fi)
figure(100);legend('Fourier Transform: real','Fourier Transform: imag')
figure(100);print -dpng fourier.png
figure(200);plot(t,p)
figure(200);legend('pulse')
figure(200);print -dpng pulse.png

