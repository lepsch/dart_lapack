void main() {      // Test program for the REAL             Level 1 CBLAS.
      // Based upon the original CBLAS test routine together with:
      // F06EAF Example Program Text
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SFAC;
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK0, CHECK1, CHECK2, CHECK3, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA             SFAC/9.765625e-4/;
      // .. Executable Statements ..
      WRITE (NOUT,99999);
      for (IC = 1; IC <= 10; IC++) { // 20
         ICASE = IC;
         header();

         // .. Initialize  PASS,  INCX,  INCY, and MODE for a new case. ..
         // .. the value 9999 for INCX, INCY or MODE will appear in the ..
         // .. detailed  output, if any, for cases  that do not involve ..
         // .. these parameters ..

         PASS = true;
         INCX = 9999;
         INCY = 9999;
         MODE = 9999;
         if (ICASE == 3) {
            check0(SFAC);
         } else if (ICASE == 7 || ICASE == 8 || ICASE == 9 || ICASE == 10) {
            check1(SFAC);
         } else if (ICASE == 1 || ICASE == 2 || ICASE == 5 || ICASE == 6) {
            check2(SFAC);
         } else if (ICASE == 4) {
            check3(SFAC);
         }
         // -- Print
         if (PASS) WRITE (NOUT,99998);
      } // 20
      STOP;

99999 FORMAT (' Real CBLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')
      }
      SUBROUTINE HEADER;
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Arrays ..
      String            L(10);
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA             L(1)/'CBLAS_SDOT '/;
      DATA             L(2)/'CBLAS_SAXPY '/;
      DATA             L(3)/'CBLAS_SROTG '/;
      DATA             L(4)/'CBLAS_SROT '/;
      DATA             L(5)/'CBLAS_SCOPY '/;
      DATA             L(6)/'CBLAS_SSWAP '/;
      DATA             L(7)/'CBLAS_SNRM2 '/;
      DATA             L(8)/'CBLAS_SASUM '/;
      DATA             L(9)/'CBLAS_SSCAL '/;
      DATA             L(10)/'CBLAS_ISAMAX'/;
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE);
      return;

99999 FORMAT (/' Test of subprogram number',I3,9X,A15)
      }
      SUBROUTINE CHECK0(SFAC);
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              SA, SB, SC, SS;
      int               K;
      // .. Local Arrays ..
      REAL              DA1(8), DATRUE(8), DB1(8), DBTRUE(8), DC1(8), DS1(8);
      // .. External Subroutines ..
      // EXTERNAL SROTGTEST, STEST1
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              DA1/0.3, 0.4, -0.3, -0.4, -0.3, 0.0, 0.0, 1.0/;
      DATA              DB1/0.4, 0.3, 0.4, 0.3, -0.4, 0.0, 1.0, 0.0/;
      DATA              DC1/0.6, 0.8, -0.6, 0.8, 0.6, 1.0, 0.0, 1.0/;
      DATA              DS1/0.8, 0.6, 0.8, -0.6, 0.8, 0.0, 1.0, 0.0/;
      DATA              DATRUE/0.5, 0.5, 0.5, -0.5, -0.5, 0.0, 1.0, 1.0/;
      DATA              DBTRUE/0.0, 0.6, 0.0, -0.6, 0.0, 0.0, 1.0, 0.0/;
      // .. Executable Statements ..

      // Compute true values which cannot be prestored
      // in decimal notation

      DBTRUE(1) = 1.0/0.6;
      DBTRUE(3) = -1.0/0.6;
      DBTRUE(5) = 1.0/0.6;

      for (K = 1; K <= 8; K++) { // 20
         // .. Set N=K for identification in output if any ..
         N = K;
         if (ICASE == 3) {
            // .. SROTGTEST ..
            if (K > 8) GO TO 40;
            SA = DA1(K);
            SB = DB1(K);
            srotgtest(SA,SB,SC,SS);
            stest1(SA,DATRUE(K),DATRUE(K),SFAC);
            stest1(SB,DBTRUE(K),DBTRUE(K),SFAC);
            stest1(SC,DC1(K),DC1(K),SFAC);
            stest1(SS,DS1(K),DS1(K),SFAC);
         } else {
            WRITE (NOUT,*) ' Shouldn''t be here in CHECK0';
            STOP;
         }
      } // 20
   40 RETURN;
      }
      SUBROUTINE CHECK1(SFAC);
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      int               I, LEN, NP1;
      // .. Local Arrays ..
      REAL              DTRUE1(5), DTRUE3(5), DTRUE5(8,5,2), DV(8,5,2), SA(10), STEMP(1), STRUE(8), SX(8);
      int               ITRUE2(5);
      // .. External Functions ..
      REAL              SASUMTEST, SNRM2TEST;
      int               ISAMAXTEST;
      // EXTERNAL SASUMTEST, SNRM2TEST, ISAMAXTEST
      // .. External Subroutines ..
      // EXTERNAL ITEST1, SSCALTEST, STEST, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              SA/0.3, -1.0, 0.0, 1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3/;
      DATA              DV/0.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.3, -0.4, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.2, -0.6, 0.3, 5.0, 5.0, 5.0, 5.0, 5.0, 0.1, -0.3, 0.5, -0.1, 6.0, 6.0, 6.0, 6.0, 0.1, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.3, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.3, 2.0, -0.4, 2.0, 2.0, 2.0, 2.0, 2.0, 0.2, 3.0, -0.6, 5.0, 0.3, 2.0, 2.0, 2.0, 0.1, 4.0, -0.3, 6.0, -0.5, 7.0, -0.1, 3.0/;
      DATA              DTRUE1/0.0, 0.3, 0.5, 0.7, 0.6/;
      DATA              DTRUE3/0.0, 0.3, 0.7, 1.1, 1.0/;
      DATA              DTRUE5/0.10, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.20, -0.60, 0.30, 5.0, 5.0, 5.0, 5.0, 5.0, 0.03, -0.09, 0.15, -0.03, 6.0, 6.0, 6.0, 6.0, 0.10, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.09, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.09, 2.0, -0.12, 2.0, 2.0, 2.0, 2.0, 2.0, 0.06, 3.0, -0.18, 5.0, 0.09, 2.0, 2.0, 2.0, 0.03, 4.0, -0.09, 6.0, -0.15, 7.0, -0.03, 3.0/;
      DATA              ITRUE2/0, 1, 2, 2, 3/;
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 80
         for (NP1 = 1; NP1 <= 5; NP1++) { // 60
            N = NP1 - 1;
            LEN = 2*MAX(N,1);
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               SX(I) = DV(I,NP1,INCX);
            } // 20

            if (ICASE == 7) {
               // .. SNRM2TEST ..
               STEMP(1) = DTRUE1(NP1);
               stest1(SNRM2TEST(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE == 8) {
               // .. SASUMTEST ..
               STEMP(1) = DTRUE3(NP1);
               stest1(SASUMTEST(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE == 9) {
               // .. SSCALTEST ..
               sscaltest(N,SA((INCX-1)*5+NP1),SX,INCX);
               for (I = 1; I <= LEN; I++) { // 40
                  STRUE(I) = DTRUE5(I,NP1,INCX);
               } // 40
               stest(LEN,SX,STRUE,STRUE,SFAC);
            } else if (ICASE == 10) {
               // .. ISAMAXTEST ..
               itest1(ISAMAXTEST(N,SX,INCX),ITRUE2(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1';
               STOP;
            }
         } // 60
      } // 80
      return;
      }
      SUBROUTINE CHECK2(SFAC);
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              SA;
      int               I, J, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      REAL              DT10X(7,4,4), DT10Y(7,4,4), DT7(4,4), DT8(7,4,4), DX1(7), DY1(7), SSIZE1(4), SSIZE2(14,2), STX(7), STY(7), SX(7), SY(7);
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      REAL              SDOTTEST;
      // EXTERNAL SDOTTEST
      // .. External Subroutines ..
      // EXTERNAL SAXPYTEST, SCOPYTEST, SSWAPTEST, STEST, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              SA/0.3/;
      DATA              INCXS/1, 2, -2, -1/;
      DATA              INCYS/1, -2, 1, -2/;
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/;
      DATA              NS/0, 1, 2, 4/;
      DATA              DX1/0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4/;
      DATA              DY1/0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8/;
      DATA              DT7/0.0, 0.30, 0.21, 0.62, 0.0, 0.30, -0.07, 0.85, 0.0, 0.30, -0.79, -0.74, 0.0, 0.30, 0.33, 1.27/;
      DATA              DT8/0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.15, 0.94, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.9, 0.48, 0.0, 0.0, 0.0, 0.0, 0.38, -0.9, 0.57, 0.7, -0.75, 0.2, 0.98, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.72, 0.0, 0.0, 0.0, 0.0, 0.0, 0.38, -0.63, 0.15, 0.88, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.7, -0.75, 0.2, 1.04/;
      DATA              DT10X/0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.8, 0.1, -0.6, 0.8, 0.3, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.7, 0.1, 0.3, 0.8, -0.9, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, -0.6, 0.8, 0.0, 0.0, 0.0/;
      DATA              DT10Y/0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.9, 0.6, 0.0, 0.0, 0.0, 0.0, -0.4, -0.9, 0.9, 0.7, -0.5, 0.2, 0.6, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.4, 0.9, -0.5, 0.6, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.7, -0.5, 0.2, 0.8/;
      DATA              SSIZE1/0.0, 0.3, 1.6, 3.2/;
      DATA              SSIZE2/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17/;
      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 120
         INCX = INCXS(KI);
         INCY = INCYS(KI);
         MX = ABS(INCX);
         MY = ABS(INCY);

         for (KN = 1; KN <= 4; KN++) { // 100
            N = NS(KN);
            KSIZE = MIN(2,KN);
            LENX = LENS(KN,MX);
            LENY = LENS(KN,MY);
            // .. Initialize all argument arrays ..
            for (I = 1; I <= 7; I++) { // 20
               SX(I) = DX1(I);
               SY(I) = DY1(I);
            } // 20

            if (ICASE == 1) {
               // .. SDOTTEST ..
               stest1(SDOTTEST(N,SX,INCX,SY,INCY),DT7(KN,KI), SSIZE1(KN),SFAC);
            } else if (ICASE == 2) {
               // .. SAXPYTEST ..
               saxpytest(N,SA,SX,INCX,SY,INCY);
               for (J = 1; J <= LENY; J++) { // 40
                  STY(J) = DT8(J,KN,KI);
               } // 40
               stest(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC);
            } else if (ICASE == 5) {
               // .. SCOPYTEST ..
               for (I = 1; I <= 7; I++) { // 60
                  STY(I) = DT10Y(I,KN,KI);
               } // 60
               scopytest(N,SX,INCX,SY,INCY);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0);
            } else if (ICASE == 6) {
               // .. SSWAPTEST ..
               sswaptest(N,SX,INCX,SY,INCY);
               for (I = 1; I <= 7; I++) { // 80
                  STX(I) = DT10X(I,KN,KI);
                  STY(I) = DT10Y(I,KN,KI);
               } // 80
               stest(LENX,SX,STX,SSIZE2(1,1),1.0);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2';
               STOP;
            }
         } // 100
      } // 120
      return;
      }
      SUBROUTINE CHECK3(SFAC);
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              SC, SS;
      int               I, K, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      REAL              COPYX(5), COPYY(5), DT9X(7,4,4), DT9Y(7,4,4), DX1(7), DY1(7), MWPC(11), MWPS(11), MWPSTX(5), MWPSTY(5), MWPTX(11,5), MWPTY(11,5), MWPX(5), MWPY(5), SSIZE2(14,2), STX(7), STY(7), SX(7), SY(7);
      int               INCXS(4), INCYS(4), LENS(4,2), MWPINX(11), MWPINY(11), MWPN(11), NS(4);
      // .. External Subroutines ..
      // EXTERNAL SROTTEST, STEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              INCXS/1, 2, -2, -1/;
      DATA              INCYS/1, -2, 1, -2/;
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/;
      DATA              NS/0, 1, 2, 4/;
      DATA              DX1/0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4/;
      DATA              DY1/0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8/;
      DATA              SC, SS/0.8, 0.6/;
      DATA              DT9X/0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, -0.22, 1.06, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.66, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.96, 0.1, -0.76, 0.8, 0.90, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.90, 0.1, -0.22, 0.8, 0.18, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, -0.76, 1.12, 0.0, 0.0, 0.0/;
      DATA              DT9Y/0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.54, 0.08, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -0.9, -0.12, 0.0, 0.0, 0.0, 0.0, 0.64, -0.9, -0.30, 0.7, -0.18, 0.2, 0.28, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -1.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.64, -1.26, 0.54, 0.20, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.7, -0.18, 0.2, 0.16/;
      DATA              SSIZE2/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17/;
      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 60
         INCX = INCXS(KI);
         INCY = INCYS(KI);
         MX = ABS(INCX);
         MY = ABS(INCY);

         for (KN = 1; KN <= 4; KN++) { // 40
            N = NS(KN);
            KSIZE = MIN(2,KN);
            LENX = LENS(KN,MX);
            LENY = LENS(KN,MY);

            if (ICASE == 4) {
               // .. SROTTEST ..
               for (I = 1; I <= 7; I++) { // 20
                  SX(I) = DX1(I);
                  SY(I) = DY1(I);
                  STX(I) = DT9X(I,KN,KI);
                  STY(I) = DT9Y(I,KN,KI);
               } // 20
               srottest(N,SX,INCX,SY,INCY,SC,SS);
               stest(LENX,SX,STX,SSIZE2(1,KSIZE),SFAC);
               stest(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK3';
               STOP;
            }
         } // 40
      } // 60

      MWPC(1) = 1;
      for (I = 2; I <= 11; I++) { // 80
         MWPC(I) = 0;
      } // 80
      MWPS(1) = 0;
      for (I = 2; I <= 6; I++) { // 100
         MWPS(I) = 1;
      } // 100
      for (I = 7; I <= 11; I++) { // 120
         MWPS(I) = -1;
      } // 120
      MWPINX(1) = 1;
      MWPINX(2) = 1;
      MWPINX(3) = 1;
      MWPINX(4) = -1;
      MWPINX(5) = 1;
      MWPINX(6) = -1;
      MWPINX(7) = 1;
      MWPINX(8) = 1;
      MWPINX(9) = -1;
      MWPINX(10) = 1;
      MWPINX(11) = -1;
      MWPINY(1) = 1;
      MWPINY(2) = 1;
      MWPINY(3) = -1;
      MWPINY(4) = -1;
      MWPINY(5) = 2;
      MWPINY(6) = 1;
      MWPINY(7) = 1;
      MWPINY(8) = -1;
      MWPINY(9) = -1;
      MWPINY(10) = 2;
      MWPINY(11) = 1;
      for (I = 1; I <= 11; I++) { // 140
         MWPN(I) = 5;
      } // 140
      MWPN(5) = 3;
      MWPN(10) = 3;
      for (I = 1; I <= 5; I++) { // 160
         MWPX(I) = I;
         MWPY(I) = I;
         MWPTX(1,I) = I;
         MWPTY(1,I) = I;
         MWPTX(2,I) = I;
         MWPTY(2,I) = -I;
         MWPTX(3,I) = 6 - I;
         MWPTY(3,I) = I - 6;
         MWPTX(4,I) = I;
         MWPTY(4,I) = -I;
         MWPTX(6,I) = 6 - I;
         MWPTY(6,I) = I - 6;
         MWPTX(7,I) = -I;
         MWPTY(7,I) = I;
         MWPTX(8,I) = I - 6;
         MWPTY(8,I) = 6 - I;
         MWPTX(9,I) = -I;
         MWPTY(9,I) = I;
         MWPTX(11,I) = I - 6;
         MWPTY(11,I) = 6 - I;
      } // 160
      MWPTX(5,1) = 1;
      MWPTX(5,2) = 3;
      MWPTX(5,3) = 5;
      MWPTX(5,4) = 4;
      MWPTX(5,5) = 5;
      MWPTY(5,1) = -1;
      MWPTY(5,2) = 2;
      MWPTY(5,3) = -2;
      MWPTY(5,4) = 4;
      MWPTY(5,5) = -3;
      MWPTX(10,1) = -1;
      MWPTX(10,2) = -3;
      MWPTX(10,3) = -5;
      MWPTX(10,4) = 4;
      MWPTX(10,5) = 5;
      MWPTY(10,1) = 1;
      MWPTY(10,2) = 2;
      MWPTY(10,3) = 2;
      MWPTY(10,4) = 4;
      MWPTY(10,5) = 3;
      for (I = 1; I <= 11; I++) { // 200
         INCX = MWPINX(I);
         INCY = MWPINY(I);
         for (K = 1; K <= 5; K++) { // 180
            COPYX(K) = MWPX(K);
            COPYY(K) = MWPY(K);
            MWPSTX(K) = MWPTX(I,K);
            MWPSTY(K) = MWPTY(I,K);
         } // 180
         srottest(MWPN(I),COPYX,INCX,COPYY,INCY,MWPC(I),MWPS(I));
         stest(5,COPYX,MWPSTX,MWPSTX,SFAC);
         stest(5,COPYY,MWPSTY,MWPSTY,SFAC);
      } // 200
      return;
      }
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC);
      // ********************************* STEST **************************

      // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
      // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
      // NEGLIGIBLE.

      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalar Arguments ..
      REAL             SFAC;
      int              LEN;
      // .. Array Arguments ..
      REAL             SCOMP(LEN), SSIZE(LEN), STRUE(LEN);
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SD;
      int              I;
      // .. External Functions ..
      REAL             SDIFF;
      // EXTERNAL SDIFF
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I);
         IF (SDIFF(ABS(SSIZE(I))+ABS(SFAC*SD),ABS(SSIZE(I))) == 0.0) GO TO 40;

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if ( !PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = false;
         WRITE (NOUT,99999);
         WRITE (NOUT,99998);
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, I, SCOMP(I), STRUE(I), SD, SSIZE(I);
      } // 40
      return;

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE  I                            ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
      }
      SUBROUTINE STEST1(SCOMP1,STRUE1,SSIZE,SFAC);
      // ************************* STEST1 *****************************

      // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
      // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
      // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      REAL              SCOMP1, SFAC, STRUE1;
      // .. Array Arguments ..
      REAL              SSIZE(*);
      // .. Local Arrays ..
      REAL              SCOMP(1), STRUE(1);
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Executable Statements ..

      SCOMP(1) = SCOMP1;
      STRUE(1) = STRUE1;
      stest(1,SCOMP,STRUE,SSIZE,SFAC);

      return;
      }
      REAL             FUNCTION SDIFF(SA,SB);
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      REAL                            SA, SB;
      // .. Executable Statements ..
      SDIFF = SA - SB;
      return;
      }
      SUBROUTINE ITEST1(ICOMP,ITRUE);
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
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      int               ID;
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..

      if (ICOMP == ITRUE) GO TO 40;

                             // HERE ICOMP IS NOT EQUAL TO ITRUE.

      if ( !PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
      PASS = false;
      WRITE (NOUT,99999);
      WRITE (NOUT,99998);
   20 ID = ICOMP - ITRUE;
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, ICOMP, ITRUE, ID;
      } // 40
      return;

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE                               ', ' COMP                                TRUE     DIFFERENCE', /1X)
99997 FORMAT (1X,I4,I3,3I5,2I36,I12)
      }
