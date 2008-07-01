#! /usr/bin/matlab

clear all

format short 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%   input: LN, Lmin, D, SD, RX, RY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
LN   = 9 ;  % # of levels in computational domain
Lmin = 0 ;  % minimum  z at level 0 
D = [  30   180   20  12.7  60  12.7   20  180  30   ] ;   % domain division
%SD= [   1    1    1   1   1    1     1   1   1       ] ;   % how many sub-elements                          
%SD= [   4    15    2    1    6    1     2   15   4   ] ;   % how many sub-elements                          
 SD= [   3    12    2    1    4    1     2   12   3   ] ;   % how many sub-elements                          
%SD= [   8    32    4    2    8    2     4   32   8  ] ;   % how many sub-elements                          
FD= [   0     1    0    1    0   -1     0   -1   0   ] ;   % taper true:1 (decrease) true:-1 (increase) false:0 
RX= [  84.76  0   50    0   30    0    50   0   84.76] ;   % radius in x
RY= [   42    0   12    0   10    0    12   0   42   ] ;   % radius in y 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%  generate files for mesh                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid1 = fopen ('levels','w'); % read file for n2to3 custumize
fid3 = fopen ('LEVEL' ,'w'); % include 'LEVEL' for usrdat2 in .usr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totl=length(D)                ;

Lsum     = Lmin               ;
fprintf(fid1,'%10.2f \n',Lmin);

ind   = 1    ;
fN(1) = Lmin ;
Dnum  = 0    ;

for i = 1:totl   
    num = SD (i)      ;
    Dnum = Dnum + D(i);
    Dsum(i) = Dnum    ; 
    for ii=1:num  
        if (i == 1) 
         dlevel= Dsum(i)/SD(i)            ;
	else
         dlevel= (Dsum(i)-Dsum(i-1))/SD(i);
	end
	tmp   =  dlevel     ;
        Lsum  =  Lsum + tmp ;
        ind   =  ind+1      ;
        fN(ind)  = Lsum     ;
        fprintf(fid1,'%10.2f \n',Lsum);
    end
end

subl= length(fN);  % total sub-levels
Lmax = fN(subl) ;  % maximum z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind   = 0     ;
Lsum  = Lmin  ;

for i = 1:totl
    num = SD (i);
    for ii=1:num  
        dlevel= D(i)/SD(i)  ;
	tmp   =  dlevel     ;
        Lsum  =  Lsum + tmp ; 
        ind   =  ind+1      ;
        fN(ind)  = Lsum     ;
    end
end

subl = length(fN)   ;
Lmax =  fN(subl)    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  declare variables for fortran subroutine 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid3,'        \n\n\n',totl);
fprintf(fid3,'        real     zl(0:%i)\n',totl);
fprintf(fid3,'        real     rx(1:%i)\n',totl);
fprintf(fid3,'        real     ry(1:%i)\n',totl);
fprintf(fid3,'        real     lmin, lmax \n');
fprintf(fid3,'        integer  totl       \n');
fprintf(fid3,'        logical  taper(1:%i)\n',totl);
fprintf(fid3,'        logical  rtaper(1:%i)\n\n\n',totl);
fprintf(fid3,'        totl =  %5i  \n',totl);
fprintf(fid3,'        lmin =  %8.2f \n',    Lmin);
fprintf(fid3,'        lmax =  %8.2f \n\n\n',Lmax);

Lsum  = Lmin;

fprintf(fid3,'        zl(0) = %10.2f\n',Lmin);

Lsum  =0    ;

for i = 1:totl
        dlevel= D(i)          ;
	tmp   =  dlevel       ;
        Lsum  =  Lsum + tmp   ;  
        fprintf(fid3,'        zl(%i) = %10.2f\n',i,Lsum);
end

        fprintf(fid3,'        \n\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  assign logical values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:totl
        dlevel= FD(i)   ;
	if (dlevel==0)  ;
        fprintf(fid3,'        TAPER(%i) = .FALSE. \n',i);
	else
        fprintf(fid3,'        TAPER(%i) = .TRUE. \n',i);
	end
end
        fprintf(fid3,'        \n\n\n');
for i = 1:totl
        dlevel= FD(i)   ;
	if (dlevel== -1); 
        fprintf(fid3,'        RTAPER(%i) = .TRUE. \n',i);
	else
        fprintf(fid3,'        RTAPER(%i) = .FALSE. \n',i);
	end
end
        fprintf(fid3,'        \n\n\n');


for i = 1:totl
        dlevel= FD(i)  ;
	if (dlevel==0) ;
	xr(i)=RX(i)/2  ;
        fprintf(fid3,'        rx(%i) = %10.2f \n',i,xr(i));
	else
	xr(i)=xr(i-1)  ;      % assign same values of previous level      
        fprintf(fid3,'        rx(%i) = %10.2f \n',i,xr(i));
	end
end
        fprintf(fid3,'        \n\n\n');
for i = 1:totl
        dlevel= FD(i)  ;
	if (dlevel==0) ;
	yr(i)=RY(i)/2  ;
        fprintf(fid3,'        ry(%i) = %10.2f \n',i,yr(i));
	else
	yr(i)=yr(i-1) ;       % assign same values of previous level      
        fprintf(fid3,'        ry(%i) = %10.2f \n',i,yr(i));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print out how many sub levels are used ?

TotolSubLevels = subl


fclose(fid1) ;
fclose(fid3) ;


