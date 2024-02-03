      void main() {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SFAC;
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK0, CHECK1, CHECK2, CHECK3, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      const SFAC = 9.765625e-4;
      // .. Executable Statements ..
      WRITE (NOUT,99999);
      for (IC = 1; IC <= 13; IC++) { // 20
         ICASE = IC;
         header();

         // .. Initialize  PASS,  INCX,  and INCY for a new case. ..
         // .. the value 9999 for INCX or INCY will appear in the ..
         // .. detailed  output, if any, for cases  that do not involve ..
         // .. these parameters ..

         PASS = true;
         INCX = 9999;
         INCY = 9999;
         if (ICASE == 3 || ICASE == 11) {
            check0(SFAC);
         } else if (ICASE == 7 || ICASE == 8 || ICASE == 9 || ICASE == 10) {
            check1(SFAC);
         } else if (ICASE == 1 || ICASE == 2 || ICASE == 5 || ICASE == 6 || ICASE == 12 || ICASE == 13) {
            check2(SFAC);
         } else if (ICASE == 4) {
            check3(SFAC);
         }
         // -- Print
         if (PASS) WRITE (NOUT,99998);
      } // 20
      STOP;

99999 FORMAT (' Real BLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')
      }
      void header() {
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
      DATA             L(1)/' SDOT '/;
      DATA             L(2)/'SAXPY '/;
      DATA             L(3)/'SROTG '/;
      DATA             L(4)/' SROT '/;
      DATA             L(5)/'SCOPY '/;
      DATA             L(6)/'SSWAP '/;
      DATA             L(7)/'SNRM2 '/;
      DATA             L(8)/'SASUM '/;
      DATA             L(9)/'SSCAL '/;
      DATA             L(10)/'ISAMAX'/;
      DATA             L(11)/'SROTMG'/;
      DATA             L(12)/'SROTM '/;
      DATA             L(13)/'SDSDOT'/;
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE);
      return;

99999 FORMAT (/' Test of subprogram number',I3,12X,A6)
      }
      void check0(SFAC) {
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              D12, SA, SB, SC, SS;
      int               I, K;
      // .. Local Arrays ..
      REAL              DA1(8), DATRUE(8), DB1(8), DBTRUE(8), DC1(8), DS1(8), DAB(4,9), DTEMP(9), DTRUE(9,9);
      // .. External Subroutines ..
      // EXTERNAL SROTG, SROTMG, STEST, STEST1
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      const DA1 = [0.3, 0.4, -0.3, -0.4, -0.3, 0.0, 0.0, 1.0];
      const DB1 = [0.4, 0.3, 0.4, 0.3, -0.4, 0.0, 1.0, 0.0];
      const DC1 = [0.6, 0.8, -0.6, 0.8, 0.6, 1.0, 0.0, 1.0];
      const DS1 = [0.8, 0.6, 0.8, -0.6, 0.8, 0.0, 1.0, 0.0];
      const DATRUE = [0.5, 0.5, 0.5, -0.5, -0.5, 0.0, 1.0, 1.0];
      const DBTRUE = [0.0, 0.6, 0.0, -0.6, 0.0, 0.0, 1.0, 0.0];
      // INPUT FOR MODIFIED GIVENS
      const DAB = [ .1,.3,1.2,.2, .7, .2, .6, 4.2, 0.0,0.0,0.0,0.0, 4.0, -1.0, 2.0, 4.0, 6e-10, 2e-2, 1e5, 10.0, 4e10, 2e-2, 1e-5, 10.0, 2e-10, 4e-2, 1e5, 10.0, 2e10, 4e-2, 1e-5, 10.0, 4.0, -2.0, 8.0, 4.0    ];
// TRUE RESULTS FOR MODIFIED GIVENS
      const DTRUE = [0.0,0.0, 1.3, .2, 0.0,0.0,0.0, .5, 0.0, 0.0,0.0, 4.5, 4.2, 1.0, .5, 0.0,0.0,0.0, 0.0,0.0,0.0,0.0, -2.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0, 4.0, -1.0, 0.0,0.0,0.0,0.0, 0.0, 15e-3, 0.0, 10.0, -1.0, 0.0, -1e-4, 0.0, 1.0, 0.0,0.0, 6144e-5, 10.0, -1.0, 4096.0, -1e6, 0.0, 1.0, 0.0,0.0,15.0,10.0,-1.0, 5e-5, 0.0,1.0,0.0, 0.0,0.0, 15.0, 10.0, -1. E0, 5e5, -4096.0, 1.0, 4096e-6, 0.0,0.0, 7.0, 4.0, 0.0,0.0, -.5, -.25, 0.0];
                    // 4096 = 2 ** 12
      const D12 = 4096.0;
      DTRUE(1,1) = 12.0 / 130.0;
      DTRUE(2,1) = 36.0 / 130.0;
      DTRUE(7,1) = -1.0 / 6.0;
      DTRUE(1,2) = 14.0 / 75.0;
      DTRUE(2,2) = 49.0 / 75.0;
      DTRUE(9,2) = 1.0 / 7.0;
      DTRUE(1,5) = 45e-11 * (D12 * D12);
      DTRUE(3,5) = 4e5 / (3.0 * D12);
      DTRUE(6,5) = 1.0 / D12;
      DTRUE(8,5) = 1e4 / (3.0 * D12);
      DTRUE(1,6) = 4e10 / (1.5 * D12 * D12);
      DTRUE(2,6) = 2e-2 / 1.5;
      DTRUE(8,6) = 5e-7 * D12;
      DTRUE(1,7) = 4.0 / 150.0;
      DTRUE(2,7) = (2e-10 / 1.5) * (D12 * D12);
      DTRUE(7,7) = -DTRUE(6,5);
      DTRUE(9,7) = 1e4 / D12;
      DTRUE(1,8) = DTRUE(1,7);
      DTRUE(2,8) = 2e10 / (1.5 * D12 * D12);
      DTRUE(1,9) = 32.0 / 7.0;
      DTRUE(2,9) = -16.0 / 7.0;
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
            // .. SROTG ..
            if (K > 8) GO TO 40;
            SA = DA1(K);
            SB = DB1(K);
            srotg(SA,SB,SC,SS);
            stest1(SA,DATRUE(K),DATRUE(K),SFAC);
            stest1(SB,DBTRUE(K),DBTRUE(K),SFAC);
            stest1(SC,DC1(K),DC1(K),SFAC);
            stest1(SS,DS1(K),DS1(K),SFAC);
         } else if (ICASE == 11) {
            // .. SROTMG ..
            for (I = 1; I <= 4; I++) {
               DTEMP(I)= DAB(I,K);
               DTEMP(I+4) = 0.0;
            }
            DTEMP(9) = 0.0;
            srotmg(DTEMP(1),DTEMP(2),DTEMP(3),DTEMP(4),DTEMP(5));
            stest(9,DTEMP,DTRUE(1,K),DTRUE(1,K),SFAC);
         } else {
            WRITE (NOUT,*) ' Shouldn''t be here in CHECK0';
            STOP;
         }
      } // 20
   40 return;
      }
      void check1(SFAC) {
      // .. Parameters ..
      int               NOUT;
      REAL              THRESH;
      const             NOUT=6, THRESH=10.0;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      int               I, IX, LEN, NP1;
      // .. Local Arrays ..
      REAL              DTRUE1(5), DTRUE3(5), DTRUE5(8,5,2), DV(8,5,2), DVR(8), SA(10), STEMP(1), STRUE(8), SX(8), SXR(15);
      int               ITRUE2(5), ITRUEC(5);
      // .. External Functions ..
      //- REAL              SASUM, SNRM2;
      //- int               ISAMAX;
      // EXTERNAL SASUM, SNRM2, ISAMAX
      // .. External Subroutines ..
      // EXTERNAL ITEST1, SB1NRM2, SSCAL, STEST, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      const SA = [0.3, -1.0, 0.0, 1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
      const DV = [0.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.3, -0.4, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.2, -0.6, 0.3, 5.0, 5.0, 5.0, 5.0, 5.0, 0.1, -0.3, 0.5, -0.1, 6.0, 6.0, 6.0, 6.0, 0.1, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.3, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.3, 2.0, -0.4, 2.0, 2.0, 2.0, 2.0, 2.0, 0.2, 3.0, -0.6, 5.0, 0.3, 2.0, 2.0, 2.0, 0.1, 4.0, -0.3, 6.0, -0.5, 7.0, -0.1, 3.0];
      const DVR = [8.0, -7.0, 9.0, 5.0, 9.0, 8.0, 7.0, 7.0];
      const DTRUE1 = [0.0, 0.3, 0.5, 0.7, 0.6];
      const DTRUE3 = [0.0, 0.3, 0.7, 1.1, 1.0];
      const DTRUE5 = [0.10, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.20, -0.60, 0.30, 5.0, 5.0, 5.0, 5.0, 5.0, 0.03, -0.09, 0.15, -0.03, 6.0, 6.0, 6.0, 6.0, 0.10, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.09, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.09, 2.0, -0.12, 2.0, 2.0, 2.0, 2.0, 2.0, 0.06, 3.0, -0.18, 5.0, 0.09, 2.0, 2.0, 2.0, 0.03, 4.0, -0.09, 6.0, -0.15, 7.0, -0.03, 3.0];
      const ITRUE2 = [0, 1, 2, 2, 3];
      const ITRUEC = [0, 1, 1, 1, 1];
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 80
         for (NP1 = 1; NP1 <= 5; NP1++) { // 60
            N = NP1 - 1;
            LEN = 2*max(N,1);
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               SX(I) = DV(I,NP1,INCX);
            } // 20

            if (ICASE == 7) {
               // .. SNRM2 ..
               // Test scaling when some entries are tiny or huge
               sb1nrm2(N,(INCX-2)*2,THRESH);
               sb1nrm2(N,INCX,THRESH);
               // Test with hardcoded mid range entries
               STEMP(1) = DTRUE1(NP1);
               stest1(SNRM2(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE == 8) {
               // .. SASUM ..
               STEMP(1) = DTRUE3(NP1);
               stest1(SASUM(N,SX,INCX),STEMP(1),STEMP,SFAC);
            } else if (ICASE == 9) {
               // .. SSCAL ..
               sscal(N,SA((INCX-1)*5+NP1),SX,INCX);
               for (I = 1; I <= LEN; I++) { // 40
                  STRUE(I) = DTRUE5(I,NP1,INCX);
               } // 40
               stest(LEN,SX,STRUE,STRUE,SFAC);
            } else if (ICASE == 10) {
               // .. ISAMAX ..
               itest1(ISAMAX(N,SX,INCX),ITRUE2(NP1));
               for (I = 1; I <= LEN; I++) { // 100
                  SX(I) = 42.0;
               } // 100
               itest1(ISAMAX(N,SX,INCX),ITRUEC(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1';
               STOP;
            }
         } // 60
         if (ICASE == 10) {
            N = 8;
            IX = 1;
            for (I = 1; I <= N; I++) { // 120
               SXR(IX) = DVR(I);
               IX = IX + INCX;
            } // 120
            itest1(ISAMAX(N,SXR,INCX),3);
         }
      } // 80
      return;
      }
      void check2(SFAC) {
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              SA;
      int               I, J, KI, KN, KNI, KPAR, KSIZE, LENX, LENY, LINCX, LINCY, MX, MY;
      // .. Local Arrays ..
      REAL              DT10X(7,4,4), DT10Y(7,4,4), DT7(4,4), DT8(7,4,4), DX1(7), DY1(7), SSIZE1(4), SSIZE2(14,2), SSIZE3(4), SSIZE(7), STX(7), STY(7), SX(7), SY(7), DPAR(5,4), DT19X(7,4,16),DT19XA(7,4,4), DT19XB(7,4,4), DT19XC(7,4,4),DT19XD(7,4,4), DT19Y(7,4,16), DT19YA(7,4,4),DT19YB(7,4,4), DT19YC(7,4,4), DT19YD(7,4,4), DTEMP(5), ST7B(4,4), STY0(1), SX0(1), SY0(1);
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      //- REAL              SDOT, SDSDOT;
      // EXTERNAL SDOT, SDSDOT
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SROTM, SSWAP, STEST, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      EQUIVALENCE (DT19X(1,1,1),DT19XA(1,1,1)),(DT19X(1,1,5), DT19XB(1,1,1)),(DT19X(1,1,9),DT19XC(1,1,1)), (DT19X(1,1,13),DT19XD(1,1,1))       EQUIVALENCE (DT19Y(1,1,1),DT19YA(1,1,1)),(DT19Y(1,1,5), DT19YB(1,1,1)),(DT19Y(1,1,9),DT19YC(1,1,1)), (DT19Y(1,1,13),DT19YD(1,1,1));

      const SA = 0.3;
      const INCXS = [1, 2, -2, -1];
      const INCYS = [1, -2, 1, -2];
      const LENS = [1, 1, 2, 4, 1, 1, 3, 7];
      const NS = [0, 1, 2, 4];
      const DX1 = [0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4];
      const DY1 = [0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8];
      const DT7 = [0.0, 0.30, 0.21, 0.62, 0.0, 0.30, -0.07, 0.85, 0.0, 0.30, -0.79, -0.74, 0.0, 0.30, 0.33, 1.27];
      const ST7B = [ .1, .4, .31, .72,     .1, .4, .03, .95, .1, .4, -.69, -.64,   .1, .4, .43, 1.37];
      const DT8 = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.15, 0.94, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.9, 0.48, 0.0, 0.0, 0.0, 0.0, 0.38, -0.9, 0.57, 0.7, -0.75, 0.2, 0.98, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.72, 0.0, 0.0, 0.0, 0.0, 0.0, 0.38, -0.63, 0.15, 0.88, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.7, -0.75, 0.2, 1.04];
      const DT10X = [0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.8, 0.1, -0.6, 0.8, 0.3, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.7, 0.1, 0.3, 0.8, -0.9, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, -0.6, 0.8, 0.0, 0.0, 0.0];
      const DT10Y = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.9, 0.6, 0.0, 0.0, 0.0, 0.0, -0.4, -0.9, 0.9, 0.7, -0.5, 0.2, 0.6, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.4, 0.9, -0.5, 0.6, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.7, -0.5, 0.2, 0.8];
      const SSIZE1 = [0.0, 0.3, 1.6, 3.2];
      const SSIZE2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17];
      const SSIZE3 = [ .1, .4, 1.7, 3.3 ];

                          // FOR DROTM

      const DPAR = [-2.0,  0.0,0.0,0.0,0.0, -1.0,  2.0, -3.0, -4.0,  5.0, 0.0,  0.0,  2.0, -3.0,  0.0, 1.0,  5.0,  2.0,  0.0, -4.0];
                         // TRUE X RESULTS F0R ROTATIONS DROTM
      const DT19XA = [.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, -.8,                  0.0,0.0,0.0,0.0,0.0,0.0, -.9,                  0.0,0.0,0.0,0.0,0.0,0.0, 3.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,   .1,             0.0,0.0,0.0,0.0,0.0, -.8,  3.8,             0.0,0.0,0.0,0.0,0.0, -.9,  2.8,             0.0,0.0,0.0,0.0,0.0, 3.5,  -.4,             0.0,0.0,0.0,0.0,0.0, .6,   .1,  -.5,   .8,          0.0,0.0,0.0, -.8,  3.8, -2.2, -1.2,          0.0,0.0,0.0, -.9,  2.8, -1.4, -1.3,          0.0,0.0,0.0, 3.5,  -.4, -2.2,  4.7,          0.0,0.0,0.0];

      const DT19XB = [.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, -.8,                  0.0,0.0,0.0,0.0,0.0,0.0, -.9,                  0.0,0.0,0.0,0.0,0.0,0.0, 3.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,   .1,  -.5,             0.0,0.0,0.0,0.0, 0.0,    .1, -3.0,             0.0,0.0,0.0,0.0, -.3,   .1, -2.0,             0.0,0.0,0.0,0.0, 3.3,   .1, -2.0,             0.0,0.0,0.0,0.0, .6,   .1,  -.5,   .8,   .9,  -.3,  -.4, -2.0,   .1,  1.4,   .8,   .6,  -.3, -2.8, -1.8,   .1,  1.3,   .8,  0.0,   -.3, -1.9, 3.8,   .1, -3.1,   .8,  4.8,  -.3, -1.5 ];

      const DT19XC = [.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, -.8,                  0.0,0.0,0.0,0.0,0.0,0.0, -.9,                  0.0,0.0,0.0,0.0,0.0,0.0, 3.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,   .1,  -.5,             0.0,0.0,0.0,0.0, 4.8,   .1, -3.0,             0.0,0.0,0.0,0.0, 3.3,   .1, -2.0,             0.0,0.0,0.0,0.0, 2.1,   .1, -2.0,             0.0,0.0,0.0,0.0, .6,   .1,  -.5,   .8,   .9,  -.3,  -.4, -1.6,   .1, -2.2,   .8,  5.4,  -.3, -2.8, -1.5,   .1, -1.4,   .8,  3.6,  -.3, -1.9, 3.7,   .1, -2.2,   .8,  3.6,  -.3, -1.5 ];

      const DT19XD = [.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,                  0.0,0.0,0.0,0.0,0.0,0.0, -.8,                  0.0,0.0,0.0,0.0,0.0,0.0, -.9,                  0.0,0.0,0.0,0.0,0.0,0.0, 3.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .6,   .1,             0.0,0.0,0.0,0.0,0.0, -.8, -1.0,             0.0,0.0,0.0,0.0,0.0, -.9,  -.8,             0.0,0.0,0.0,0.0,0.0, 3.5,   .8,             0.0,0.0,0.0,0.0,0.0, .6,   .1,  -.5,   .8,          0.0,0.0,0.0, -.8, -1.0,  1.4, -1.6,          0.0,0.0,0.0, -.9,  -.8,  1.3, -1.6,          0.0,0.0,0.0, 3.5,   .8, -3.1,  4.8,          0.0,0.0,0.0];
                         // TRUE Y RESULTS FOR ROTATIONS DROTM
      const DT19YA = [.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .7,                  0.0,0.0,0.0,0.0,0.0,0.0, 1.7,                  0.0,0.0,0.0,0.0,0.0,0.0, -2.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,  -.9,             0.0,0.0,0.0,0.0,0.0, .7, -4.8,             0.0,0.0,0.0,0.0,0.0, 1.7,  -.7,             0.0,0.0,0.0,0.0,0.0, -2.6,  3.5,             0.0,0.0,0.0,0.0,0.0, .5,  -.9,   .3,   .7,          0.0,0.0,0.0, .7, -4.8,  3.0,  1.1,          0.0,0.0,0.0, 1.7,  -.7,  -.7,  2.3,          0.0,0.0,0.0, -2.6,  3.5,  -.7, -3.6,          0.0,0.0,0.0];

      const DT19YB = [.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .7,                  0.0,0.0,0.0,0.0,0.0,0.0, 1.7,                  0.0,0.0,0.0,0.0,0.0,0.0, -2.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,  -.9,   .3,             0.0,0.0,0.0,0.0, 4.0,  -.9,  -.3,             0.0,0.0,0.0,0.0, -.5,  -.9,  1.5,             0.0,0.0,0.0,0.0, -1.5,  -.9, -1.8,             0.0,0.0,0.0,0.0, .5,  -.9,   .3,   .7,  -.6,   .2,   .8, 3.7,  -.9, -1.2,   .7, -1.5,   .2,  2.2, -.3,  -.9,  2.1,   .7, -1.6,   .2,  2.0, -1.6,  -.9, -2.1,   .7,  2.9,   .2, -3.8 ];

      const DT19YC = [.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .7,                  0.0,0.0,0.0,0.0,0.0,0.0, 1.7,                  0.0,0.0,0.0,0.0,0.0,0.0, -2.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,  -.9,             0.0,0.0,0.0,0.0,0.0, 4.0, -6.3,             0.0,0.0,0.0,0.0,0.0, -.5,   .3,             0.0,0.0,0.0,0.0,0.0, -1.5,  3.0,             0.0,0.0,0.0,0.0,0.0, .5,  -.9,   .3,   .7,          0.0,0.0,0.0, 3.7, -7.2,  3.0,  1.7,          0.0,0.0,0.0, -.3,   .9,  -.7,  1.9,          0.0,0.0,0.0, -1.6,  2.7,  -.7, -3.4,          0.0,0.0,0.0];

      const DT19YD = [.5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,                  0.0,0.0,0.0,0.0,0.0,0.0, .7,                  0.0,0.0,0.0,0.0,0.0,0.0, 1.7,                  0.0,0.0,0.0,0.0,0.0,0.0, -2.6,                  0.0,0.0,0.0,0.0,0.0,0.0, .5,  -.9,   .3,             0.0,0.0,0.0,0.0, .7,  -.9,  1.2,             0.0,0.0,0.0,0.0, 1.7,  -.9,   .5,             0.0,0.0,0.0,0.0, -2.6,  -.9, -1.3,             0.0,0.0,0.0,0.0, .5,  -.9,   .3,   .7,  -.6,   .2,   .8, .7,  -.9,  1.2,   .7, -1.5,   .2,  1.6, 1.7,  -.9,   .5,   .7, -1.6,   .2,  2.4, -2.6,  -.9, -1.3,   .7,  2.9,   .2, -4.0 ];

      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 120
         INCX = INCXS(KI);
         INCY = INCYS(KI);
         MX = (INCX).abs();
         MY = (INCY).abs();

         for (KN = 1; KN <= 4; KN++) { // 100
            N = NS(KN);
            KSIZE = min(2,KN);
            LENX = LENS(KN,MX);
            LENY = LENS(KN,MY);
            // .. Initialize all argument arrays ..
            for (I = 1; I <= 7; I++) { // 20
               SX(I) = DX1(I);
               SY(I) = DY1(I);
            } // 20

            if (ICASE == 1) {
               // .. SDOT ..
               stest1(SDOT(N,SX,INCX,SY,INCY),DT7(KN,KI),SSIZE1(KN) ,SFAC);
            } else if (ICASE == 2) {
               // .. SAXPY ..
               saxpy(N,SA,SX,INCX,SY,INCY);
               for (J = 1; J <= LENY; J++) { // 40
                  STY(J) = DT8(J,KN,KI);
               } // 40
               stest(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC);
            } else if (ICASE == 5) {
               // .. SCOPY ..
               for (I = 1; I <= 7; I++) { // 60
                  STY(I) = DT10Y(I,KN,KI);
               } // 60
               scopy(N,SX,INCX,SY,INCY);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0);
               if (KI == 1) {
                  SX0(1) = 42.0;
                  SY0(1) = 43.0;
                  if (N == 0) {
                     STY0(1) = SY0(1);
                  } else {
                     STY0(1) = SX0(1);
                  }
                  LINCX = INCX;
                  INCX = 0;
                  LINCY = INCY;
                  INCY = 0;
                  scopy(N,SX0,INCX,SY0,INCY);
                  stest(1,SY0,STY0,SSIZE2(1,1),1.0);
                  INCX = LINCX;
                  INCY = LINCY;
               }
            } else if (ICASE == 6) {
               // .. SSWAP ..
               sswap(N,SX,INCX,SY,INCY);
               for (I = 1; I <= 7; I++) { // 80
                  STX(I) = DT10X(I,KN,KI);
                  STY(I) = DT10Y(I,KN,KI);
               } // 80
               stest(LENX,SX,STX,SSIZE2(1,1),1.0);
               stest(LENY,SY,STY,SSIZE2(1,1),1.0);
            } else if (ICASE == 12) {
               // .. SROTM ..
               KNI=KN+4*(KI-1);
               for (KPAR = 1; KPAR <= 4; KPAR++) {
                  for (I = 1; I <= 7; I++) {
                     SX(I) = DX1(I);
                     SY(I) = DY1(I);
                     STX(I)= DT19X(I,KPAR,KNI);
                     STY(I)= DT19Y(I,KPAR,KNI);
                  }

                  for (I = 1; I <= 5; I++) {
                     DTEMP(I) = DPAR(I,KPAR);
                  }

                  for (I = 1; I <= LENX; I++) {
                     SSIZE(I)=STX(I);
                  }
                    // SEE REMARK ABOVE ABOUT DT11X(1,2,7)
                        // AND DT11X(5,3,8).
                  if ((KPAR == 2) && (KNI == 7)) SSIZE(1) = 2.4;
                  IF ((KPAR == 3) && (KNI == 8)) SSIZE(5) = 1.8;

                  srotm(N,SX,INCX,SY,INCY,DTEMP);
                  stest(LENX,SX,STX,SSIZE,SFAC);
                  stest(LENY,SY,STY,STY,SFAC);
               }
            } else if (ICASE == 13) {
               // .. SDSROT ..
               stest1(SDSDOT(N,.1,SX,INCX,SY,INCY), ST7B(KN,KI),SSIZE3(KN),SFAC);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2';
               STOP;
            }
         } // 100
      } // 120
      return;
      }
      void check3(SFAC) {
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, N;
      bool              PASS;
      // .. Local Scalars ..
      REAL              SC, SS;
      int               I, K, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      REAL              COPYX(5), COPYY(5), DT9X(7,4,4), DT9Y(7,4,4), DX1(7), DY1(7), MWPC(11), MWPS(11), MWPSTX(5), MWPSTY(5), MWPTX(11,5), MWPTY(11,5), MWPX(5), MWPY(5), SSIZE2(14,2), STX(7), STY(7), SX(7), SY(7);
      int               INCXS(4), INCYS(4), LENS(4,2), MWPINX(11), MWPINY(11), MWPN(11), NS(4);
      // .. External Subroutines ..
      // EXTERNAL SROT, STEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Data statements ..
      const INCXS = [1, 2, -2, -1];
      const INCYS = [1, -2, 1, -2];
      const LENS = [1, 1, 2, 4, 1, 1, 3, 7];
      const NS = [0, 1, 2, 4];
      const DX1 = [0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4];
      const DY1 = [0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8];
      final (SC, SS) = (0.8, 0.6);
      const DT9X = [0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, -0.22, 1.06, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.66, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.96, 0.1, -0.76, 0.8, 0.90, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.90, 0.1, -0.22, 0.8, 0.18, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, -0.76, 1.12, 0.0, 0.0, 0.0];
      const DT9Y = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.54, 0.08, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -0.9, -0.12, 0.0, 0.0, 0.0, 0.0, 0.64, -0.9, -0.30, 0.7, -0.18, 0.2, 0.28, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -1.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.64, -1.26, 0.54, 0.20, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.7, -0.18, 0.2, 0.16];
      const SSIZE2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17];
      // .. Executable Statements ..

      for (KI = 1; KI <= 4; KI++) { // 60
         INCX = INCXS(KI);
         INCY = INCYS(KI);
         MX = (INCX).abs();
         MY = (INCY).abs();

         for (KN = 1; KN <= 4; KN++) { // 40
            N = NS(KN);
            KSIZE = min(2,KN);
            LENX = LENS(KN,MX);
            LENY = LENS(KN,MY);

            if (ICASE == 4) {
               // .. SROT ..
               for (I = 1; I <= 7; I++) { // 20
                  SX(I) = DX1(I);
                  SY(I) = DY1(I);
                  STX(I) = DT9X(I,KN,KI);
                  STY(I) = DT9Y(I,KN,KI);
               } // 20
               srot(N,SX,INCX,SY,INCY,SC,SS);
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
         srot(MWPN(I),COPYX,INCX,COPYY,INCY,MWPC(I),MWPS(I));
         stest(5,COPYX,MWPSTX,MWPSTX,SFAC);
         stest(5,COPYY,MWPSTY,MWPSTY,SFAC);
      } // 200
      return;
      }
      void stest(LEN,SCOMP,STRUE,SSIZE,SFAC) {
      // ********************************* STEST **************************

      // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
      // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
      // NEGLIGIBLE.

      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int              NOUT;
      REAL             ZERO;
      const            NOUT=6, ZERO=0.0;
      // .. Scalar Arguments ..
      REAL             SFAC;
      int              LEN;
      // .. Array Arguments ..
      REAL             SCOMP(LEN), SSIZE(LEN), STRUE(LEN);
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SD;
      int              I;
      // .. External Functions ..
      //- REAL             SDIFF;
      // EXTERNAL SDIFF
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I);
         if ((SFAC*SD).abs() <= (SSIZE(I)).abs()*EPSILON(ZERO)) GO TO 40;

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if ( !PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = false;
         WRITE (NOUT,99999);
         WRITE (NOUT,99998);
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, I, SCOMP(I), STRUE(I), SD, SSIZE(I);
      } // 40
      return;

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY  I                            ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,2I5,I3,2E36.8,2E12.4)
      }
      void stest1(SCOMP1,STRUE1,SSIZE,SFAC) {
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
      REAL sdiff(SA,SB) {
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      REAL                            SA, SB;
      // .. Executable Statements ..
      SDIFF = SA - SB;
      return;
      }
      void itest1(ICOMP,ITRUE) {
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

      if (ICOMP == ITRUE) GO TO 40;

                             // HERE ICOMP IS NOT EQUAL TO ITRUE.

      if ( !PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
      PASS = false;
      WRITE (NOUT,99999);
      WRITE (NOUT,99998);
   20 ID = ICOMP - ITRUE;
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, ICOMP, ITRUE, ID;
      } // 40
      return;

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY                               ', ' COMP                                TRUE     DIFFERENCE', /1X)
99997 FORMAT (1X,I4,I3,2I5,2I36,I12)
      }
      void sb1nrm2(N,INCX,THRESH) {
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

      // IMPLICIT NONE
      // .. Scalar Arguments ..
      int               INCX, N;
      REAL              THRESH;

// =====================================================================
      // .. Parameters ..
      int               NMAX, NOUT, NV;
      const             NMAX=20, NOUT=6, NV=10;
      REAL              HALF, ONE, TWO, ZERO;
      const             HALF=0.5, ONE=1.0, TWO= 2.0, ZERO=0.0;
      // .. External Functions ..
      //- REAL              SNRM2;
      // EXTERNAL SNRM2
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // .. Model parameters ..
      REAL              BIGNUM, SAFMAX, SAFMIN, SMLNUM, ULP;
      const             BIGNUM=0.1014120480e+32, SAFMAX=0.8507059173e+38, SAFMIN=0.1175494351e-37, SMLNUM=0.9860761315e-31, ULP=0.1192092896e-06;
      // .. Local Scalars ..
      REAL              ROGUE, SNRM, TRAT, V0, V1, WORKSSQ, Y1, Y2, YMAX, YMIN, YNRM, ZNRM;
      int               I, IV, IW, IX;
      bool              FIRST;
      // .. Local Arrays ..
      REAL              VALUES(NV), WORK(NMAX), X(NMAX), Z(NMAX);
      // .. Executable Statements ..
      VALUES(1) = ZERO;
      VALUES(2) = TWO*SAFMIN;
      VALUES(3) = SMLNUM;
      VALUES(4) = ULP;
      VALUES(5) = ONE;
      VALUES(6) = ONE / ULP;
      VALUES(7) = BIGNUM;
      VALUES(8) = SAFMAX;
      VALUES(9) = SXVALS(V0,2);
      VALUES(10) = SXVALS(V0,3);
      ROGUE = -1234.5678;
      FIRST = true;

      // Check that the arrays are large enough

      if (N*(INCX).abs() > NMAX) {
         WRITE (NOUT,99) "SNRM2", NMAX, INCX, N, N*(INCX).abs();
         return;
      }

      // Zero-sized inputs are tested in STEST1.
      if (N <= 0) {
         return;
      }

      // Generate (N-1) values in (-1,1).

      for (I = 2; I <= N; I++) {
         random_number(WORK(I));
         WORK(I) = ONE - TWO*WORK(I);
      }

      // Compute the sum of squares of the random values
      // by an unscaled algorithm.

      WORKSSQ = ZERO;
      for (I = 2; I <= N; I++) {
         WORKSSQ = WORKSSQ + WORK(I)*WORK(I);
      }

      // Construct the test vector with one known value
      // and the rest from the random work array multiplied
      // by a scaling factor.

      for (IV = 1; IV <= NV; IV++) {
         V0 = VALUES(IV);
         if ((V0).abs() > ONE) {
         V0 = V0*HALF;
         }
         Z(1) = V0;
         for (IW = 1; IW <= NV; IW++) {
            V1 = VALUES(IW);
            if ((V1).abs() > ONE) {
               V1 = (V1*HALF) / sqrt(REAL(N));
            }
            for (I = 2; I <= N; I++) {
               Z(I) = V1*WORK(I);
            }

            // Compute the expected value of the 2-norm

            Y1 = (V0).abs();
            if (N > 1) {
               Y2 = (V1).abs()*sqrt(WORKSSQ);
            } else {
               Y2 = ZERO;
            }
            YMIN = min(Y1, Y2);
            YMAX = max(Y1, Y2);

            // Expected value is NaN if either is NaN. The test
            // for YMIN == YMAX avoids further computation if both
            // are infinity.

            if ((Y1 != Y1) || (Y2 != Y2)) {
               // add to propagate NaN
               YNRM = Y1 + Y2;
            } else if (YMIN == YMAX) {
               YNRM = sqrt(TWO)*YMAX;
            } else if (YMAX == ZERO) {
               YNRM = ZERO;
            } else {
               YNRM = YMAX*sqrt(ONE + (YMIN / YMAX)**2);
            }

            // Fill the input array to SNRM2 with steps of incx

            for (I = 1; I <= N; I++) {
               X(I) = ROGUE;
            }
            IX = 1;
            if (INCX < 0) IX = 1 - (N-1)*INCX;
            for (I = 1; I <= N; I++) {
               X(IX) = Z(I);
               IX = IX + INCX;
            }

            // Call SNRM2 to compute the 2-norm

            SNRM = SNRM2(N,X,INCX);

            // Compare SNRM and ZNRM.  Roundoff error grows like O(n)
            // in this implementation so we scale the test ratio accordingly.

            if (INCX == 0) {
               ZNRM = sqrt(REAL(N))*(X(1)).abs();
            } else {
               ZNRM = YNRM;
            }

            // The tests for NaN rely on the compiler not being overly
            // aggressive and removing the statements altogether.
            if ((SNRM != SNRM) || (ZNRM != ZNRM)) {
               if ((SNRM != SNRM).NEQV.(ZNRM != ZNRM)) {
                  TRAT = ONE / ULP;
               } else {
                  TRAT = ZERO;
               }
            } else if (SNRM == ZNRM) {
               TRAT = ZERO;
            } else if (ZNRM == ZERO) {
               TRAT = SNRM / ULP;
            } else {
               TRAT = ((SNRM-ZNRM).abs() / ZNRM) / (REAL(N)*ULP);
            }
            if ((TRAT != TRAT) || (TRAT >= THRESH)) {
               if (FIRST) {
                  FIRST = false;
                  WRITE(NOUT,99999);
               }
               WRITE (NOUT,98) "SNRM2", N, INCX, IV, IW, TRAT;
            }
         }
      }
99999 FORMAT ('                                       FAIL')
   99 FORMAT ( ' Not enough space to test ', A6, ': NMAX = ',I6, ', INCX = ',I6,/,'   N = ',I6,', must be at least ',I6 );
   98 FORMAT( 1X, A6, ': N=', I6,', INCX=', I4, ', IV=', I2, ', IW=', I2, ', test=', E15.8 );
      return;
      CONTAINS;
      REAL sxvals(XX,K) {
      // .. Scalar Arguments ..
      REAL              XX;
      int               K;
      // .. Local Scalars ..
      REAL              X, Y, YY, Z;
      // .. Intrinsic Functions ..
      // INTRINSIC HUGE
      // .. Executable Statements ..
      Y = HUGE(XX);
      Z = YY;
      if (K == 1) {
         X = -Z;
      } else if (K == 2) {
         X = Z;
      } else if (K == 3) {
         X = Z / Z;
      }
      SXVALS = X;
      return;
      }
      }
