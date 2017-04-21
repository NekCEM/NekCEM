      SUBROUTINE CONVERT(C,T)
C
      REAL*4 T
      CHARACTER*4  C
C
      CHARACTER*12 VC
      CHARACTER*1  V1(12)
      EQUIVALENCE (V1,VC)
      CHARACTER*1  CSIGN
C
      REAL*4 W
      CHARACTER*4  WC
      CHARACTER*1  W1(4)
      EQUIVALENCE (W1,W)
      EQUIVALENCE (W1,WC)
C
      CHARACTER*1 ALPH64(0:63)
      SAVE        ALPH64
      DATA        ALPH64 
     $  / '1','2','3','4','5','6','7','8','9','0'
     $   ,'a','b','c','d','e','f','g','h','i','j'
     $   ,'k','l','m','n','o','p','q','r','s','t'
     $   ,'u','v','w','x','y','z'
     $   ,'A','B','C','D','E','F','G','H','I','J'
     $   ,'K','L','M','N','O','P','Q','R','S','T'
     $   ,'U','V','W','X','Y','Z','+','-'/
C
C     Find out the usual decimal format for T
C
      WRITE(VC,10) T
   10 FORMAT(E12.5)
C
C     Begin converting the mantissa to base 64
C
      READ(VC,11) MANTIS
   11 FORMAT(3X,I5)
C
C     Exponent
C
      READ(VC,21) IEXP
   21 FORMAT(9X,I3)
      IF (IEXP.LT.-30) THEN
         IEXP=0
         MANTIS=0
      ELSEIF (IEXP.GT.30) THEN
         IEXP=25
      ENDIF
C     We assume that the exponent is bounded by 31.
      IEXP=IEXP+31
C
C     Sign?
C
      READ(VC,12) CSIGN
   12 FORMAT(A1)
      IF (CSIGN.EQ.'-') MANTIS=-MANTIS
      MANTIS=MANTIS+131072
C
C     ONES,TENS, HUNDREDS
C
      IONE=MOD(MANTIS,64)
      ITMP=MANTIS/64
      ITEN=MOD(ITMP,64)
      ITMP=ITMP/64
      IHUN=MOD(ITMP,64)
C
C     Compute alpha equivalent
C
      W1(1)=ALPH64(IHUN)
      W1(2)=ALPH64(ITEN)
      W1(3)=ALPH64(IONE)
      W1(4)=ALPH64(IEXP)
C
C     Convert the input value
      C=WC
      RETURN
      END
