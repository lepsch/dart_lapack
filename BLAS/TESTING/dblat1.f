void main() {
*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Scalars ..
      double           SFAC;
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK0, CHECK1, CHECK2, CHECK3, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      DATA             SFAC/9.765625D-4/
      // .. Executable Statements ..
      WRITE (NOUT,99999)
      for (IC = 1; IC <= 13; IC++) { // 20
         ICASE = IC
         header();

         // .. Initialize  PASS,  INCX,  and INCY for a new case. ..
         // .. the value 9999 for INCX or INCY will appear in the ..
         // .. detailed  output, if any, for cases  that do not involve ..
         // .. these parameters ..

         PASS = .TRUE.
         INCX = 9999
         INCY = 9999
         if (ICASE.EQ.3 .OR. ICASE.EQ.11) {
            check0(SFAC);
         } else if (ICASE.EQ.7 .OR. ICASE.EQ.8 .OR. ICASE.EQ.9 .OR. ICASE.EQ.10) {
            check1(SFAC);
         } else if (ICASE.EQ.1 .OR. ICASE.EQ.2 .OR. ICASE.EQ.5 .OR. ICASE.EQ.6 .OR. ICASE.EQ.12 .OR. ICASE.EQ.13) {
            check2(SFAC);
         } else if (ICASE.EQ.4) {
            check3(SFAC);
         }
         // -- Print
         if (PASS) WRITE (NOUT,99998);
      } // 20
      STOP

99999 FORMAT (' Real BLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')

      // End of DBLAT1

      }
      SUBROUTINE HEADER
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Arrays ..
      String           L(13);
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      DATA             L(1)/' DDOT '/
      DATA             L(2)/'DAXPY '/
      DATA             L(3)/'DROTG '/
      DATA             L(4)/' DROT '/
      DATA             L(5)/'DCOPY '/
      DATA             L(6)/'DSWAP '/
      DATA             L(7)/'DNRM2 '/
      DATA             L(8)/'DASUM '/
      DATA             L(9)/'DSCAL '/
      DATA             L(10)/'IDAMAX'/
      DATA             L(11)/'DROTMG'/
      DATA             L(12)/'DROTM '/
      DATA             L(13)/'DSDOT '/
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE)
      RETURN

99999 FORMAT (/' Test of subprogram number',I3,12X,A6)

      // End of HEADER

      }
      SUBROUTINE CHECK0(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      double            SA, SB, SC, SS, D12;
      int               I, K;
      // .. Local Arrays ..
      double            DA1(8), DATRUE(8), DB1(8), DBTRUE(8), DC1(8), DS1(8), DAB(4,9), DTEMP(9), DTRUE(9,9);
      // .. External Subroutines ..
      // EXTERNAL DROTG, DROTMG, STEST, STEST1
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      DATA              DA1/0.3D0, 0.4D0, -0.3D0, -0.4D0, -0.3D0, 0.0D0, 0.0D0, 1.0D0/
      DATA              DB1/0.4D0, 0.3D0, 0.4D0, 0.3D0, -0.4D0, 0.0D0, 1.0D0, 0.0D0/
      DATA              DC1/0.6D0, 0.8D0, -0.6D0, 0.8D0, 0.6D0, 1.0D0, 0.0D0, 1.0D0/
      DATA              DS1/0.8D0, 0.6D0, 0.8D0, -0.6D0, 0.8D0, 0.0D0, 1.0D0, 0.0D0/
      DATA              DATRUE/0.5D0, 0.5D0, 0.5D0, -0.5D0, -0.5D0, 0.0D0, 1.0D0, 1.0D0/
      DATA              DBTRUE/0.0D0, 0.6D0, 0.0D0, -0.6D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0/
      // INPUT FOR MODIFIED GIVENS
      DATA DAB/ .1D0,.3D0,1.2D0,.2D0, .7D0, .2D0, .6D0, 4.2D0, 0.D0,0.D0,0.D0,0.D0, 4.D0, -1.D0, 2.D0, 4.D0, 6.D-10, 2.D-2, 1.D5, 10.D0, 4.D10, 2.D-2, 1.D-5, 10.D0, 2.D-10, 4.D-2, 1.D5, 10.D0, 2.D10, 4.D-2, 1.D-5, 10.D0, 4.D0, -2.D0, 8.D0, 4.D0    /
*    TRUE RESULTS FOR MODIFIED GIVENS
      DATA DTRUE/0.D0,0.D0, 1.3D0, .2D0, 0.D0,0.D0,0.D0, .5D0, 0.D0, 0.D0,0.D0, 4.5D0, 4.2D0, 1.D0, .5D0, 0.D0,0.D0,0.D0, 0.D0,0.D0,0.D0,0.D0, -2.D0, 0.D0,0.D0,0.D0,0.D0, 0.D0,0.D0,0.D0, 4.D0, -1.D0, 0.D0,0.D0,0.D0,0.D0, 0.D0, 15.D-3, 0.D0, 10.D0, -1.D0, 0.D0, -1.D-4, 0.D0, 1.D0, 0.D0,0.D0, 6144.D-5, 10.D0, -1.D0, 4096.D0, -1.D6, 0.D0, 1.D0, 0.D0,0.D0,15.D0,10.D0,-1.D0, 5.D-5, 0.D0,1.D0,0.D0, 0.D0,0.D0, 15.D0, 10.D0, -1. D0, 5.D5, -4096.D0, 1.D0, 4096.D-6, 0.D0,0.D0, 7.D0, 4.D0, 0.D0,0.D0, -.5D0, -.25D0, 0.D0/
                    // 4096 = 2 ** 12
      DATA D12  /4096.D0/
      DTRUE(1,1) = 12.D0 / 130.D0
      DTRUE(2,1) = 36.D0 / 130.D0
      DTRUE(7,1) = -1.D0 / 6.D0
      DTRUE(1,2) = 14.D0 / 75.D0
      DTRUE(2,2) = 49.D0 / 75.D0
      DTRUE(9,2) = 1.D0 / 7.D0
      DTRUE(1,5) = 45.D-11 * (D12 * D12)
      DTRUE(3,5) = 4.D5 / (3.D0 * D12)
      DTRUE(6,5) = 1.D0 / D12
      DTRUE(8,5) = 1.D4 / (3.D0 * D12)
      DTRUE(1,6) = 4.D10 / (1.5D0 * D12 * D12)
      DTRUE(2,6) = 2.D-2 / 1.5D0
      DTRUE(8,6) = 5.D-7 * D12
      DTRUE(1,7) = 4.D0 / 150.D0
      DTRUE(2,7) = (2.D-10 / 1.5D0) * (D12 * D12)
      DTRUE(7,7) = -DTRUE(6,5)
      DTRUE(9,7) = 1.D4 / D12
      DTRUE(1,8) = DTRUE(1,7)
      DTRUE(2,8) = 2.D10 / (1.5D0 * D12 * D12)
      DTRUE(1,9) = 32.D0 / 7.D0
      DTRUE(2,9) = -16.D0 / 7.D0
      // .. Executable Statements ..

      // Compute true values which cannot be prestored
      // in decimal notation

      DBTRUE(1) = 1.0D0/0.6D0
      DBTRUE(3) = -1.0D0/0.6D0
      DBTRUE(5) = 1.0D0/0.6D0

      for (K = 1; K <= 8; K++) { // 20
         // .. Set N=K for identification in output if any ..
         N = K
         if (ICASE.EQ.3) {
            // .. DROTG ..
            if (K.GT.8) GO TO 40;
            SA = DA1(K)
            SB = DB1(K)
            drotg(SA,SB,SC,SS);
            stest1(SA,DATRUE(K),DATRUE(K),SFAC);
            stest1(SB,DBTRUE(K),DBTRUE(K),SFAC);
            stest1(SC,DC1(K),DC1(K),SFAC);
            stest1(SS,DS1(K),DS1(K),SFAC);
         } else if (ICASE.EQ.11) {
            // .. DROTMG ..
            for (I = 1; I <= 4; I++) {
               DTEMP(I)= DAB(I,K)
               DTEMP(I+4) = 0.0
            }
            DTEMP(9) = 0.0
            drotmg(DTEMP(1),DTEMP(2),DTEMP(3),DTEMP(4),DTEMP(5));
            stest(9,DTEMP,DTRUE(1,K),DTRUE(1,K),SFAC);
         } else {
            WRITE (NOUT,*) ' Shouldn''t be here in CHECK0'
            STOP
         }
      } // 20
   40 RETURN

      // End of CHECK0

      }
      SUBROUTINE CHECK1(SFAC)
      // .. Parameters ..
      double            THRESH;
      int               NOUT;
      const             NOUT=6, THRESH=10.0D0;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      int               I, IX, LEN, NP1;
      // .. Local Arrays ..
      double            DTRUE1(5), DTRUE3(5), DTRUE5(8,5,2), DV(8,5,2), DVR(8), SA(10), STEMP(1), STRUE(8), SX(8), SXR(15);
      int               ITRUE2(5), ITRUEC(5);
      // .. External Functions ..
      double            DASUM, DNRM2;
      int               IDAMAX;
      // EXTERNAL DASUM, DNRM2, IDAMAX
      // .. External Subroutines ..
      // EXTERNAL ITEST1, DB1NRM2, DSCAL, STEST, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      DATA              SA/0.3D0, -1.0D0, 0.0D0, 1.0D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0/
      DATA              DV/0.1D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 0.3D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 0.3D0, -0.4D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 0.2D0, -0.6D0, 0.3D0, 5.0D0, 5.0D0, 5.0D0, 5.0D0, 5.0D0, 0.1D0, -0.3D0, 0.5D0, -0.1D0, 6.0D0, 6.0D0, 6.0D0, 6.0D0, 0.1D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 0.3D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 0.3D0, 2.0D0, -0.4D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 0.2D0, 3.0D0, -0.6D0, 5.0D0, 0.3D0, 2.0D0, 2.0D0, 2.0D0, 0.1D0, 4.0D0, -0.3D0, 6.0D0, -0.5D0, 7.0D0, -0.1D0, 3.0D0/
      DATA              DVR/8.0D0, -7.0D0, 9.0D0, 5.0D0, 9.0D0, 8.0D0, 7.0D0, 7.0D0/
      DATA              DTRUE1/0.0D0, 0.3D0, 0.5D0, 0.7D0, 0.6D0/
      DATA              DTRUE3/0.0D0, 0.3D0, 0.7D0, 1.1D0, 1.0D0/
      DATA              DTRUE5/0.10D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, -0.3D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 3.0D0, 0.0D0, 0.0D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 4.0D0, 0.20D0, -0.60D0, 0.30D0, 5.0D0, 5.0D0, 5.0D0, 5.0D0, 5.0D0, 0.03D0, -0.09D0, 0.15D0, -0.03D0, 6.0D0, 6.0D0, 6.0D0, 6.0D0, 0.10D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 8.0D0, 0.09D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 9.0D0, 0.09D0, 2.0D0, -0.12D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 2.0D0, 0.06D0, 3.0D0, -0.18D0, 5.0D0, 0.09D0, 2.0D0, 2.0D0, 2.0D0, 0.03D0, 4.0D0, -0.09D0, 6.0D0, -0.15D0, 7.0D0, -0.03D0, 3.0D0/
      DATA              ITRUE2/0, 1, 2, 2, 3/
      DATA              ITRUEC/0, 1, 1, 1, 1/
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 80
         for (NP1 = 1; NP1 <= 5; NP1++) { // 60
            N = NP1 - 1
            LEN = 2*MAX(N,1)
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               SX(I) = DV(I,NP1,INCX)
            } // 20

            if (ICASE.EQ.7) {
               // .. DNRM2 ..
               // Test scaling when some entries are tiny or huge
               db1nrm2(N,(INCX-2)*2,THRESH);
               db1nrm2(N,INCX,THRESH);
               // Test with hardcoded mid range entries
               STEMP(1) = DTRUE1(NP1)
               stest1(DNRM2(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE.EQ.8) {
               // .. DASUM ..
               STEMP(1) = DTRUE3(NP1)
               stest1(DASUM(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE.EQ.9) {
               // .. DSCAL ..
               dscal(N,SA((INCX-1)*5+NP1),SX,INCX);
               for (I = 1; I <= LEN; I++) { // 40
                  STRUE(I) = DTRUE5(I,NP1,INCX)
               } // 40
               stest(LEN,SX,STRUE,STRUE,SFAC);
            } else if (ICASE.EQ.10) {
               // .. IDAMAX ..
               itest1(IDAMAX(N,SX,INCX),ITRUE2(NP1));
               for (I = 1; I <= LEN; I++) { // 100
                  SX(I) = 42.0D0
               } // 100
               itest1(IDAMAX(N,SX,INCX),ITRUEC(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            }
         } // 60
         if (ICASE.EQ.10) {
            N = 8
            IX = 1
            for (I = 1; I <= N; I++) { // 120
               SXR(IX) = DVR(I)
               IX = IX + INCX
            } // 120
            itest1(IDAMAX(N,SXR,INCX),3);
         }
      } // 80
      RETURN

      // End of CHECK1

      }
      SUBROUTINE CHECK2(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      double            SA;
      int               I, J, KI, KN, KNI, KPAR, KSIZE, LENX, LENY, LINCX, LINCY, MX, MY;
      // .. Local Arrays ..
      double            DT10X(7,4,4), DT10Y(7,4,4), DT7(4,4), DT8(7,4,4), DX1(7), DY1(7), SSIZE1(4), SSIZE2(14,2), SSIZE(7), STX(7), STY(7), SX(7), SY(7), DPAR(5,4), DT19X(7,4,16),DT19XA(7,4,4), DT19XB(7,4,4), DT19XC(7,4,4),DT19XD(7,4,4), DT19Y(7,4,16), DT19YA(7,4,4),DT19YB(7,4,4), DT19YC(7,4,4), DT19YD(7,4,4), DTEMP(5), STY0(1), SX0(1), SY0(1);
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      double            DDOT, DSDOT;
      // EXTERNAL DDOT, DSDOT
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DROTM, DSWAP, STEST, STEST1, TESTDSDOT
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      EQUIVALENCE (DT19X(1,1,1),DT19XA(1,1,1)),(DT19X(1,1,5), DT19XB(1,1,1)),(DT19X(1,1,9),DT19XC(1,1,1)), (DT19X(1,1,13),DT19XD(1,1,1))       EQUIVALENCE (DT19Y(1,1,1),DT19YA(1,1,1)),(DT19Y(1,1,5), DT19YB(1,1,1)),(DT19Y(1,1,9),DT19YC(1,1,1)), (DT19Y(1,1,13),DT19YD(1,1,1))

      DATA              SA/0.3D0/
      DATA              INCXS/1, 2, -2, -1/
      DATA              INCYS/1, -2, 1, -2/
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/
      DATA              NS/0, 1, 2, 4/
      DATA              DX1/0.6D0, 0.1D0, -0.5D0, 0.8D0, 0.9D0, -0.3D0, -0.4D0/
      DATA              DY1/0.5D0, -0.9D0, 0.3D0, 0.7D0, -0.6D0, 0.2D0, 0.8D0/
      DATA              DT7/0.0D0, 0.30D0, 0.21D0, 0.62D0, 0.0D0, 0.30D0, -0.07D0, 0.85D0, 0.0D0, 0.30D0, -0.79D0, -0.74D0, 0.0D0, 0.30D0, 0.33D0, 1.27D0/
      DATA              DT8/0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, -0.87D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, -0.87D0, 0.15D0, 0.94D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.35D0, -0.9D0, 0.48D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.38D0, -0.9D0, 0.57D0, 0.7D0, -0.75D0, 0.2D0, 0.98D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.35D0, -0.72D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.38D0, -0.63D0, 0.15D0, 0.88D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, -0.9D0, 0.33D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.68D0, -0.9D0, 0.33D0, 0.7D0, -0.75D0, 0.2D0, 1.04D0/
      DATA              DT10X/0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, -0.9D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, -0.9D0, 0.3D0, 0.7D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.3D0, 0.1D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.8D0, 0.1D0, -0.6D0, 0.8D0, 0.3D0, -0.3D0, 0.5D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.9D0, 0.1D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.7D0, 0.1D0, 0.3D0, 0.8D0, -0.9D0, -0.3D0, 0.5D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.3D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.3D0, -0.6D0, 0.8D0, 0.0D0, 0.0D0, 0.0D0/
      DATA              DT10Y/0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.1D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.1D0, -0.5D0, 0.8D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.5D0, -0.9D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.4D0, -0.9D0, 0.9D0, 0.7D0, -0.5D0, 0.2D0, 0.6D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.5D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.4D0, 0.9D0, -0.5D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, -0.9D0, 0.1D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, -0.9D0, 0.1D0, 0.7D0, -0.5D0, 0.2D0, 0.8D0/
      DATA              SSIZE1/0.0D0, 0.3D0, 1.6D0, 3.2D0/
      DATA              SSIZE2/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0/

                          // FOR DROTM

      DATA DPAR/-2.D0,  0.D0,0.D0,0.D0,0.D0, -1.D0,  2.D0, -3.D0, -4.D0,  5.D0, 0.D0,  0.D0,  2.D0, -3.D0,  0.D0, 1.D0,  5.D0,  2.D0,  0.D0, -4.D0/
                         // TRUE X RESULTS F0R ROTATIONS DROTM
      DATA DT19XA/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0,  3.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,  2.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,  -.4D0,             0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0, -.8D0,  3.8D0, -2.2D0, -1.2D0,          0.D0,0.D0,0.D0, -.9D0,  2.8D0, -1.4D0, -1.3D0,          0.D0,0.D0,0.D0, 3.5D0,  -.4D0, -2.2D0,  4.7D0,          0.D0,0.D0,0.D0/

      DATA DT19XB/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0, 0.D0,    .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0, -.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0, 3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0, -2.0D0,   .1D0,  1.4D0,   .8D0,   .6D0,  -.3D0, -2.8D0, -1.8D0,   .1D0,  1.3D0,   .8D0,  0.D0,   -.3D0, -1.9D0, 3.8D0,   .1D0, -3.1D0,   .8D0,  4.8D0,  -.3D0, -1.5D0 /

      DATA DT19XC/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0, 4.8D0,   .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0, 3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0, 2.1D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0, -1.6D0,   .1D0, -2.2D0,   .8D0,  5.4D0,  -.3D0, -2.8D0, -1.5D0,   .1D0, -1.4D0,   .8D0,  3.6D0,  -.3D0, -1.9D0, 3.7D0,   .1D0, -2.2D0,   .8D0,  3.6D0,  -.3D0, -1.5D0 /

      DATA DT19XD/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -.8D0, -1.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -.9D0,  -.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0, 3.5D0,   .8D0,             0.D0,0.D0,0.D0,0.D0,0.D0, .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0, -.8D0, -1.0D0,  1.4D0, -1.6D0,          0.D0,0.D0,0.D0, -.9D0,  -.8D0,  1.3D0, -1.6D0,          0.D0,0.D0,0.D0, 3.5D0,   .8D0, -3.1D0,  4.8D0,          0.D0,0.D0,0.D0/
                         // TRUE Y RESULTS FOR ROTATIONS DROTM
      DATA DT19YA/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0, .7D0, -4.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0, 1.7D0,  -.7D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -2.6D0,  3.5D0,             0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0, .7D0, -4.8D0,  3.0D0,  1.1D0,          0.D0,0.D0,0.D0, 1.7D0,  -.7D0,  -.7D0,  2.3D0,          0.D0,0.D0,0.D0, -2.6D0,  3.5D0,  -.7D0, -3.6D0,          0.D0,0.D0,0.D0/

      DATA DT19YB/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0, 4.0D0,  -.9D0,  -.3D0,             0.D0,0.D0,0.D0,0.D0, -.5D0,  -.9D0,  1.5D0,             0.D0,0.D0,0.D0,0.D0, -1.5D0,  -.9D0, -1.8D0,             0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0, 3.7D0,  -.9D0, -1.2D0,   .7D0, -1.5D0,   .2D0,  2.2D0, -.3D0,  -.9D0,  2.1D0,   .7D0, -1.6D0,   .2D0,  2.0D0, -1.6D0,  -.9D0, -2.1D0,   .7D0,  2.9D0,   .2D0, -3.8D0 /

      DATA DT19YC/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0, 4.0D0, -6.3D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -.5D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,0.D0, -1.5D0,  3.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0, 3.7D0, -7.2D0,  3.0D0,  1.7D0,          0.D0,0.D0,0.D0, -.3D0,   .9D0,  -.7D0,  1.9D0,          0.D0,0.D0,0.D0, -1.6D0,  2.7D0,  -.7D0, -3.4D0,          0.D0,0.D0,0.D0/

      DATA DT19YD/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, 1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0, .7D0,  -.9D0,  1.2D0,             0.D0,0.D0,0.D0,0.D0, 1.7D0,  -.9D0,   .5D0,             0.D0,0.D0,0.D0,0.D0, -2.6D0,  -.9D0, -1.3D0,             0.D0,0.D0,0.D0,0.D0, .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0, .7D0,  -.9D0,  1.2D0,   .7D0, -1.5D0,   .2D0,  1.6D0, 1.7D0,  -.9D0,   .5D0,   .7D0, -1.6D0,   .2D0,  2.4D0, -2.6D0,  -.9D0, -1.3D0,   .7D0,  2.9D0,   .2D0, -4.0D0 /

      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 120
         INCX = INCXS(KI)
         INCY = INCYS(KI)
         MX = ABS(INCX)
         MY = ABS(INCY)

         for (KN = 1; KN <= 4; KN++) { // 100
            N = NS(KN)
            KSIZE = MIN(2,KN)
            LENX = LENS(KN,MX)
            LENY = LENS(KN,MY)
            // .. Initialize all argument arrays ..
            for (I = 1; I <= 7; I++) { // 20
               SX(I) = DX1(I)
               SY(I) = DY1(I)
            } // 20

            if (ICASE.EQ.1) {
               // .. DDOT ..
               stest1(DDOT(N,SX,INCX,SY,INCY),DT7(KN,KI),SSIZE1(KN) ,SFAC);
            } else if (ICASE.EQ.2) {
               // .. DAXPY ..
               daxpy(N,SA,SX,INCX,SY,INCY);
               for (J = 1; J <= LENY; J++) { // 40
                  STY(J) = DT8(J,KN,KI)
               } // 40
               stest(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC);
            } else if (ICASE.EQ.5) {
               // .. DCOPY ..
               for (I = 1; I <= 7; I++) { // 60
                  STY(I) = DT10Y(I,KN,KI)
               } // 60
               dcopy(N,SX,INCX,SY,INCY);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0D0);
               if (KI.EQ.1) {
                  SX0(1) = 42.0D0
                  SY0(1) = 43.0D0
                  if (N.EQ.0) {
                     STY0(1) = SY0(1)
                  } else {
                     STY0(1) = SX0(1)
                  }
                  LINCX = INCX
                  INCX = 0
                  LINCY = INCY
                  INCY = 0
                  dcopy(N,SX0,INCX,SY0,INCY);
                  stest(1,SY0,STY0,SSIZE2(1,1),1.0D0);
                  INCX = LINCX
                  INCY = LINCY
               }
            } else if (ICASE.EQ.6) {
               // .. DSWAP ..
               dswap(N,SX,INCX,SY,INCY);
               for (I = 1; I <= 7; I++) { // 80
                  STX(I) = DT10X(I,KN,KI)
                  STY(I) = DT10Y(I,KN,KI)
               } // 80
               stest(LENX,SX,STX,SSIZE2(1,1),1.0D0);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0D0);
            } else if (ICASE.EQ.12) {
               // .. DROTM ..
               KNI=KN+4*(KI-1)
               for (KPAR = 1; KPAR <= 4; KPAR++) {
                  for (I = 1; I <= 7; I++) {
                     SX(I) = DX1(I)
                     SY(I) = DY1(I)
                     STX(I)= DT19X(I,KPAR,KNI)
                     STY(I)= DT19Y(I,KPAR,KNI)
                  }

                  for (I = 1; I <= 5; I++) {
                     DTEMP(I) = DPAR(I,KPAR)
                  }

                  for (I = 1; I <= LENX; I++) {
                     SSIZE(I)=STX(I)
                  }
                    // SEE REMARK ABOVE ABOUT DT11X(1,2,7)
                        // AND DT11X(5,3,8).
                  IF ((KPAR .EQ. 2) .AND. (KNI .EQ. 7)) SSIZE(1) = 2.4D0                   IF ((KPAR .EQ. 3) .AND. (KNI .EQ. 8)) SSIZE(5) = 1.8D0

                  drotm(N,SX,INCX,SY,INCY,DTEMP);
                  stest(LENX,SX,STX,SSIZE,SFAC);
                  stest(LENY,SY,STY,STY,SFAC);
               }
            } else if (ICASE.EQ.13) {
               // .. DSDOT ..
            testdsdot(REAL(DSDOT(N,REAL(SX),INCX,REAL(SY),INCY)), REAL(DT7(KN,KI)),REAL(SSIZE1(KN)), .3125E-1);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            }
         } // 100
      } // 120
      RETURN

      // End of CHECK2

      }
      SUBROUTINE CHECK3(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      double            SC, SS;
      int               I, K, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      double            COPYX(5), COPYY(5), DT9X(7,4,4), DT9Y(7,4,4), DX1(7), DY1(7), MWPC(11), MWPS(11), MWPSTX(5), MWPSTY(5), MWPTX(11,5), MWPTY(11,5), MWPX(5), MWPY(5), SSIZE2(14,2), STX(7), STY(7), SX(7), SY(7);
      int               INCXS(4), INCYS(4), LENS(4,2), MWPINX(11), MWPINY(11), MWPN(11), NS(4);
      // .. External Subroutines ..
      // EXTERNAL DROT, STEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      DATA              INCXS/1, 2, -2, -1/
      DATA              INCYS/1, -2, 1, -2/
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/
      DATA              NS/0, 1, 2, 4/
      DATA              DX1/0.6D0, 0.1D0, -0.5D0, 0.8D0, 0.9D0, -0.3D0, -0.4D0/
      DATA              DY1/0.5D0, -0.9D0, 0.3D0, 0.7D0, -0.6D0, 0.2D0, 0.8D0/
      DATA              SC, SS/0.8D0, 0.6D0/
      DATA              DT9X/0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, -0.46D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, -0.46D0, -0.22D0, 1.06D0, 0.0D0, 0.0D0, 0.0D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.66D0, 0.1D0, -0.1D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.96D0, 0.1D0, -0.76D0, 0.8D0, 0.90D0, -0.3D0, -0.02D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -0.06D0, 0.1D0, -0.1D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.90D0, 0.1D0, -0.22D0, 0.8D0, 0.18D0, -0.3D0, -0.02D0, 0.6D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.26D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.78D0, 0.26D0, -0.76D0, 1.12D0, 0.0D0, 0.0D0, 0.0D0/
      DATA              DT9Y/0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, -0.78D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, -0.78D0, 0.54D0, 0.08D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.7D0, -0.9D0, -0.12D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.64D0, -0.9D0, -0.30D0, 0.7D0, -0.18D0, 0.2D0, 0.28D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.7D0, -1.08D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.64D0, -1.26D0, 0.54D0, 0.20D0, 0.0D0, 0.0D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, -0.9D0, 0.18D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.04D0, -0.9D0, 0.18D0, 0.7D0, -0.18D0, 0.2D0, 0.16D0/
      DATA              SSIZE2/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0/
      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 60
         INCX = INCXS(KI)
         INCY = INCYS(KI)
         MX = ABS(INCX)
         MY = ABS(INCY)

         for (KN = 1; KN <= 4; KN++) { // 40
            N = NS(KN)
            KSIZE = MIN(2,KN)
            LENX = LENS(KN,MX)
            LENY = LENS(KN,MY)

            if (ICASE.EQ.4) {
               // .. DROT ..
               for (I = 1; I <= 7; I++) { // 20
                  SX(I) = DX1(I)
                  SY(I) = DY1(I)
                  STX(I) = DT9X(I,KN,KI)
                  STY(I) = DT9Y(I,KN,KI)
               } // 20
               drot(N,SX,INCX,SY,INCY,SC,SS);
               stest(LENX,SX,STX,SSIZE2(1,KSIZE),SFAC);
               stest(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK3'
               STOP
            }
         } // 40
      } // 60

      MWPC(1) = 1
      for (I = 2; I <= 11; I++) { // 80
         MWPC(I) = 0
      } // 80
      MWPS(1) = 0
      for (I = 2; I <= 6; I++) { // 100
         MWPS(I) = 1
      } // 100
      for (I = 7; I <= 11; I++) { // 120
         MWPS(I) = -1
      } // 120
      MWPINX(1) = 1
      MWPINX(2) = 1
      MWPINX(3) = 1
      MWPINX(4) = -1
      MWPINX(5) = 1
      MWPINX(6) = -1
      MWPINX(7) = 1
      MWPINX(8) = 1
      MWPINX(9) = -1
      MWPINX(10) = 1
      MWPINX(11) = -1
      MWPINY(1) = 1
      MWPINY(2) = 1
      MWPINY(3) = -1
      MWPINY(4) = -1
      MWPINY(5) = 2
      MWPINY(6) = 1
      MWPINY(7) = 1
      MWPINY(8) = -1
      MWPINY(9) = -1
      MWPINY(10) = 2
      MWPINY(11) = 1
      for (I = 1; I <= 11; I++) { // 140
         MWPN(I) = 5
      } // 140
      MWPN(5) = 3
      MWPN(10) = 3
      for (I = 1; I <= 5; I++) { // 160
         MWPX(I) = I
         MWPY(I) = I
         MWPTX(1,I) = I
         MWPTY(1,I) = I
         MWPTX(2,I) = I
         MWPTY(2,I) = -I
         MWPTX(3,I) = 6 - I
         MWPTY(3,I) = I - 6
         MWPTX(4,I) = I
         MWPTY(4,I) = -I
         MWPTX(6,I) = 6 - I
         MWPTY(6,I) = I - 6
         MWPTX(7,I) = -I
         MWPTY(7,I) = I
         MWPTX(8,I) = I - 6
         MWPTY(8,I) = 6 - I
         MWPTX(9,I) = -I
         MWPTY(9,I) = I
         MWPTX(11,I) = I - 6
         MWPTY(11,I) = 6 - I
      } // 160
      MWPTX(5,1) = 1
      MWPTX(5,2) = 3
      MWPTX(5,3) = 5
      MWPTX(5,4) = 4
      MWPTX(5,5) = 5
      MWPTY(5,1) = -1
      MWPTY(5,2) = 2
      MWPTY(5,3) = -2
      MWPTY(5,4) = 4
      MWPTY(5,5) = -3
      MWPTX(10,1) = -1
      MWPTX(10,2) = -3
      MWPTX(10,3) = -5
      MWPTX(10,4) = 4
      MWPTX(10,5) = 5
      MWPTY(10,1) = 1
      MWPTY(10,2) = 2
      MWPTY(10,3) = 2
      MWPTY(10,4) = 4
      MWPTY(10,5) = 3
      for (I = 1; I <= 11; I++) { // 200
         INCX = MWPINX(I)
         INCY = MWPINY(I)
         for (K = 1; K <= 5; K++) { // 180
            COPYX(K) = MWPX(K)
            COPYY(K) = MWPY(K)
            MWPSTX(K) = MWPTX(I,K)
            MWPSTY(K) = MWPTY(I,K)
         } // 180
         drot(MWPN(I),COPYX,INCX,COPYY,INCY,MWPC(I),MWPS(I));
         stest(5,COPYX,MWPSTX,MWPSTX,SFAC);
         stest(5,COPYY,MWPSTY,MWPSTY,SFAC);
      } // 200
      RETURN

      // End of CHECK3

      }
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC)
      // ********************************* STEST **************************

      // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
      // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
      // NEGLIGIBLE.

      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int              NOUT;
      double           ZERO;
      const            NOUT=6, ZERO=0.0D0;
      // .. Scalar Arguments ..
      double           SFAC;
      int              LEN;
      // .. Array Arguments ..
      double           SCOMP(LEN), SSIZE(LEN), STRUE(LEN);
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Scalars ..
      double           SD;
      int              I;
      // .. External Functions ..
      double           SDIFF;
      // EXTERNAL SDIFF
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I)
         IF (ABS(SFAC*SD) .LE. ABS(SSIZE(I))*EPSILON(ZERO)) GO TO 40

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = .FALSE.
         WRITE (NOUT,99999)
         WRITE (NOUT,99998)
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, I, SCOMP(I), STRUE(I), SD, SSIZE(I)
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY  I                            ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,2I5,I3,2D36.8,2D12.4)

      // End of STEST

      }
      SUBROUTINE TESTDSDOT(SCOMP,STRUE,SSIZE,SFAC)
      // ********************************* STEST **************************

      // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
      // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
      // NEGLIGIBLE.

      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int              NOUT;
      REAL             ZERO
      const            NOUT=6, ZERO=0.0E0;
      // .. Scalar Arguments ..
      REAL             SFAC, SCOMP, SSIZE, STRUE
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SD
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Executable Statements ..

         SD = SCOMP - STRUE
         IF (ABS(SFAC*SD) .LE. ABS(SSIZE) * EPSILON(ZERO)) GO TO 40

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = .FALSE.
         WRITE (NOUT,99999)
         WRITE (NOUT,99998)
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, SCOMP, STRUE, SD, SSIZE
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY                           ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,1I5,I3,2E36.8,2E12.4)

      // End of TESTDSDOT

      }
      SUBROUTINE STEST1(SCOMP1,STRUE1,SSIZE,SFAC)
      // ************************* STEST1 *****************************

      // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
      // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
      // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      double            SCOMP1, SFAC, STRUE1;
      // .. Array Arguments ..
      double            SSIZE(*);
      // .. Local Arrays ..
      double            SCOMP(1), STRUE(1);
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Executable Statements ..

      SCOMP(1) = SCOMP1
      STRUE(1) = STRUE1
      stest(1,SCOMP,STRUE,SSIZE,SFAC);

      RETURN

      // End of STEST1

      }
      double           FUNCTION SDIFF(SA,SB);
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      double                          SA, SB;
      // .. Executable Statements ..
      SDIFF = SA - SB
      RETURN

      // End of SDIFF

      }
      SUBROUTINE ITEST1(ICOMP,ITRUE)
      // ********************************* ITEST1 *************************

      // THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
      // EQUALITY.
      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      int               ICOMP, ITRUE;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      int               ID;
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Executable Statements ..

      if (ICOMP.EQ.ITRUE) GO TO 40;

                             // HERE ICOMP IS NOT EQUAL TO ITRUE.

      if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
      PASS = .FALSE.
      WRITE (NOUT,99999)
      WRITE (NOUT,99998)
   20 ID = ICOMP - ITRUE
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, ICOMP, ITRUE, ID
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY                               ', ' COMP                                TRUE     DIFFERENCE', /1X)
99997 FORMAT (1X,I4,I3,2I5,2I36,I12)

      // End of ITEST1

      }
      SUBROUTINE DB1NRM2(N,INCX,THRESH)
      // Compare NRM2 with a reference computation using combinations
      // of the following values:

      // 0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN

      // one of these values is used to initialize x(1) and x(2:N) is
      // filled with random values from [-1,1] scaled by another of
      // these values.

      // This routine is adapted from the test suite provided by
      // Anderson E. (2017)
      // Algorithm 978: Safe Scaling in the Level 1 BLAS
      // ACM Trans Math Softw 44:1--28
      // https://doi.org/10.1145/3061665

      // .. Scalar Arguments ..
      int               INCX, N;
      double            THRESH;

*  =====================================================================
      // .. Parameters ..
      int               NMAX, NOUT, NV;
      const             NMAX=20, NOUT=6, NV=10;
      double            HALF, ONE, TWO, ZERO;
      const             HALF=0.5D+0, ONE=1.0D+0, TWO= 2.0D+0, ZERO=0.0D+0;
      // .. External Functions ..
      double            DNRM2;
      // EXTERNAL DNRM2
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // .. Model parameters ..
      double            BIGNUM, SAFMAX, SAFMIN, SMLNUM, ULP;
      const             BIGNUM=0.99792015476735990583D+292, SAFMAX=0.44942328371557897693D+308, SAFMIN=0.22250738585072013831D-307, SMLNUM=0.10020841800044863890D-291, ULP=0.22204460492503130808D-015;
      // .. Local Scalars ..
      double            ROGUE, SNRM, TRAT, V0, V1, WORKSSQ, Y1, Y2, YMAX, YMIN, YNRM, ZNRM;
      int               I, IV, IW, IX;
      bool              FIRST;
      // .. Local Arrays ..
      double            VALUES(NV), WORK(NMAX), X(NMAX), Z(NMAX);
      // .. Executable Statements ..
      VALUES(1) = ZERO
      VALUES(2) = TWO*SAFMIN
      VALUES(3) = SMLNUM
      VALUES(4) = ULP
      VALUES(5) = ONE
      VALUES(6) = ONE / ULP
      VALUES(7) = BIGNUM
      VALUES(8) = SAFMAX
      VALUES(9) = DXVALS(V0,2)
      VALUES(10) = DXVALS(V0,3)
      ROGUE = -1234.5678D+0
      FIRST = .TRUE.

      // Check that the arrays are large enough

      if (N*ABS(INCX).GT.NMAX) {
         WRITE (NOUT,99) "DNRM2", NMAX, INCX, N, N*ABS(INCX)
         RETURN
      }

      // Zero-sized inputs are tested in STEST1.
      if (N.LE.0) {
         RETURN
      }

      // Generate (N-1) values in (-1,1).

      for (I = 2; I <= N; I++) {
         random_number(WORK(I));
         WORK(I) = ONE - TWO*WORK(I)
      }

      // Compute the sum of squares of the random values
      // by an unscaled algorithm.

      WORKSSQ = ZERO
      for (I = 2; I <= N; I++) {
         WORKSSQ = WORKSSQ + WORK(I)*WORK(I)
      }

      // Construct the test vector with one known value
      // and the rest from the random work array multiplied
      // by a scaling factor.

      for (IV = 1; IV <= NV; IV++) {
         V0 = VALUES(IV)
         if (ABS(V0).GT.ONE) {
            V0 = V0*HALF
         }
         Z(1) = V0
         for (IW = 1; IW <= NV; IW++) {
            V1 = VALUES(IW)
            if (ABS(V1).GT.ONE) {
               V1 = (V1*HALF) / SQRT(DBLE(N))
            }
            for (I = 2; I <= N; I++) {
               Z(I) = V1*WORK(I)
            }

            // Compute the expected value of the 2-norm

            Y1 = ABS(V0)
            if (N.GT.1) {
               Y2 = ABS(V1)*SQRT(WORKSSQ)
            } else {
               Y2 = ZERO
            }
            YMIN = MIN(Y1, Y2)
            YMAX = MAX(Y1, Y2)

            // Expected value is NaN if either is NaN. The test
            // for YMIN == YMAX avoids further computation if both
            // are infinity.

            if ((Y1.NE.Y1).OR.(Y2.NE.Y2)) {
               // Add to propagate NaN
               YNRM = Y1 + Y2
            } else if (YMAX == ZERO) {
               YNRM = ZERO
            } else if (YMIN == YMAX) {
               YNRM = SQRT(TWO)*YMAX
            } else {
               YNRM = YMAX*SQRT(ONE + (YMIN / YMAX)**2)
            }

            // Fill the input array to DNRM2 with steps of incx

            for (I = 1; I <= N; I++) {
               X(I) = ROGUE
            }
            IX = 1
            if (INCX.LT.0) IX = 1 - (N-1)*INCX;
            for (I = 1; I <= N; I++) {
               X(IX) = Z(I)
               IX = IX + INCX
            }

            // Call DNRM2 to compute the 2-norm

            SNRM = DNRM2(N,X,INCX)

            // Compare SNRM and ZNRM.  Roundoff error grows like O(n)
            // in this implementation so we scale the test ratio accordingly.

            if (INCX.EQ.0) {
               ZNRM = SQRT(DBLE(N))*ABS(X(1))
            } else {
               ZNRM = YNRM
            }

            // The tests for NaN rely on the compiler not being overly
            // aggressive and removing the statements altogether.
            if ((SNRM.NE.SNRM).OR.(ZNRM.NE.ZNRM)) {
               if ((SNRM.NE.SNRM).NEQV.(ZNRM.NE.ZNRM)) {
                  TRAT = ONE / ULP
               } else {
                  TRAT = ZERO
               }
            } else if (SNRM == ZNRM) {
               TRAT = ZERO
            } else if (ZNRM == ZERO) {
               TRAT = SNRM / ULP
            } else {
               TRAT = (ABS(SNRM-ZNRM) / ZNRM) / (DBLE(N)*ULP)
            }
            if ((TRAT.NE.TRAT).OR.(TRAT.GE.THRESH)) {
               if (FIRST) {
                  FIRST = .FALSE.
                  WRITE(NOUT,99999)
               }
               WRITE (NOUT,98) "DNRM2", N, INCX, IV, IW, TRAT
            }
         }
      }
99999 FORMAT ('                                       FAIL')
   99 FORMAT ( ' Not enough space to test ', A6, ': NMAX = ',I6, ', INCX = ',I6,/,'   N = ',I6,', must be at least ',I6 )
   98 FORMAT( 1X, A6, ': N=', I6,', INCX=', I4, ', IV=', I2, ', IW=', I2, ', test=', E15.8 )
      RETURN
      CONTAINS
      double           FUNCTION DXVALS(XX,K);
      // .. Scalar Arguments ..
      double            XX;
      int               K;
      // .. Local Scalars ..
      double            X, Y, YY, Z;
      // .. Intrinsic Functions ..
      // INTRINSIC HUGE
      // .. Executable Statements ..
      Y = HUGE(XX)
      Z = YY
      if (K.EQ.1) {
         X = -Z
      } else if (K.EQ.2) {
         X = Z
      } else if (K.EQ.3) {
         X = Z / Z
      }
      DXVALS = X
      RETURN
      }
      }
