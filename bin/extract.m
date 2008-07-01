clear all

file1 = input('old rea file=');
file2 = input('new rea file=');
elementnum = input('total element number=');
direction = input('direction:: top(1) = 0 or bottom(-1) =1 ? ');

input_unit  = fopen (char(file1),'r') 
output_unit = fopen (char(file2),'w');

eltnum = elementnum;
preline= 129  ;
eltcon = 7    ;
cntelt = 0    ;
totelt = 0    ;
endline= eltnum*eltcon+preline;
count  = 0    ;
zline  = []   ;
elttag = char('            ELEMENT'); 

if     (direction == 0)
    zvalue = char('   1.0000000       1.0000000       1.0000000       1.000000')
elseif (direction == 1)
    zvalue = char('  -1.0000000      -1.0000000      -1.0000000      -1.000000')            
else
    stop
end

for    i = 1:endline           
       line = fgets(input_unit) ; 
       if  (i> preline) 
          lineno = i-preline    ;
	  if ((mod(lineno,7)==4)|(mod(lineno,7)==0))
	      zcheck = line(1:59);
	      if (zcheck==zvalue)
	       zline = [zline lineno]; 
	       totelt = totelt +1 ;
              end  
          end
       end
end

fclose(input_unit)

input_unit  = fopen (char(file1),'r');

for    i = 1:endline           

       line    = fgets(input_unit) ; 

       count   = count +1;

       if     ((count>=1)&(count <= 2))
              fprintf(output_unit,line) 
       elseif (count==3)
	      fprintf(output_unit,'2 DIMENSIONAL RUN\n') 
       elseif ((count>=4)&(count <= 127))
              fprintf(output_unit,line) 
       elseif (count==128)
	      fprintf(output_unit,'**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.\n') 
       elseif (count==129)
	      fprintf(output_unit,'  %6d    2  %6d            NEL,NDIM,NELV\n',totelt,totelt) 
       else

          cntglo= count-129;

          min01 = min(abs(cntglo-(zline-3)));
          min02 = min(abs(cntglo-(zline-6)));
          min11 = min(abs(cntglo-(zline-1)));
          min12 = min(abs(cntglo-(zline-2)));

          eltcheck=line(1:19);

	  if (min01==0|(min02)==0)
	     if (mod(cntglo,7)==1)
	         cntelt=cntelt+1;
                 fprintf(output_unit,'            ELEMENT  %4d [    1 ]    GROUP     0\n',cntelt) 
	     end
          end
          if (min11==0|min12==0)
               fprintf(output_unit,line) 
          end  

       end 
end

fclose(input_unit)
