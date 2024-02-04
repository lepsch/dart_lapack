      void main() {
      // Test program for the COMPLEX    Level 1 CBLAS.
      // Based upon the original CBLAS test routine together with:
      // F06GAF Example Program Text
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      double             SFAC;
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK1, CHECK2, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      const SFAC = 9.765625e-4;
      // .. Executable Statements ..
      WRITE (NOUT,99999);
      for (IC = 1; IC <= 10; IC++) { // 20
         ICASE = IC;
         header();

         // Initialize PASS, INCX, INCY, and MODE for a new case.
         // The value 9999 for INCX, INCY or MODE will appear in the
         // detailed  output, if any, for cases that do not involve
         // these parameters.

         PASS = true;
         INCX = 9999;
         INCY = 9999;
         MODE = 9999;
         if (ICASE <= 5) {
            check2(SFAC);
         } else if (ICASE >= 6) {
            check1(SFAC);
         }
         // -- Print
         if (PASS) WRITE (NOUT,99998);
      } // 20
      STOP;

99999 FORMAT (' Complex CBLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')
      }
      void header() {
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
      DATA             L(1)/'CBLAS_CDOTC'/;
      DATA             L(2)/'CBLAS_CDOTU'/;
      DATA             L(3)/'CBLAS_CAXPY'/;
      DATA             L(4)/'CBLAS_CCOPY'/;
      DATA             L(5)/'CBLAS_CSWAP'/;
      DATA             L(6)/'CBLAS_SCNRM2'/;
      DATA             L(7)/'CBLAS_SCASUM'/;
      DATA             L(8)/'CBLAS_CSCAL'/;
      DATA             L(9)/'CBLAS_CSSCAL'/;
      DATA             L(10)/'CBLAS_ICAMAX'/;
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE);
      return;

99999 FORMAT (/' Test of subprogram number',I3,9X,A15)
      }
      void check1(SFAC) {
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      Complex           CA;
      double              SA;
      int               I, J, LEN, NP1;
      // .. Local Arrays ..
      Complex           CTRUE5(8,5,2), CTRUE6(8,5,2), CV(8,5,2), CX(8), MWPCS(5), MWPCT(5);
      double              STRUE2(5), STRUE4(5);
      int               ITRUE3(5);
      // .. External Functions ..
      //- REAL              SCASUMTEST, SCNRM2TEST;
      //- int               ICAMAXTEST;
      // EXTERNAL SCASUMTEST, SCNRM2TEST, ICAMAXTEST
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSSCALTEST, CTEST, ITEST1, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      final (SA, CA) = (0.3, (0.4,-0.7));
      DATA              ((CV(I,J,1),I=1,8),J=1,5)/(0.1,0.1), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (0.3,-0.4), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (0.1,-0.3), (0.5,-0.1), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (0.1,0.1), (-0.6,0.1), (0.1,-0.3), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (0.3,0.1), (0.1,0.4), (0.4,0.1), (0.1,0.2), (2.0,3.0), (2.0,3.0), (2.0,3.0), (2.0,3.0)/;
      DATA              ((CV(I,J,2),I=1,8),J=1,5)/(0.1,0.1), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (0.3,-0.4), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (0.1,-0.3), (8.0,9.0), (0.5,-0.1), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (0.1,0.1), (3.0,6.0), (-0.6,0.1), (4.0,7.0), (0.1,-0.3), (7.0,2.0), (7.0,2.0), (7.0,2.0), (0.3,0.1), (5.0,8.0), (0.1,0.4), (6.0,9.0), (0.4,0.1), (8.0,3.0), (0.1,0.2), (9.0,4.0)/;
      const STRUE2 = [0.0, 0.5, 0.6, 0.7, 0.7];
      const STRUE4 = [0.0, 0.7, 1.0, 1.3, 1.7];
      DATA              ((CTRUE5(I,J,1),I=1,8),J=1,5)/(0.1,0.1), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (-0.16,-0.37), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (-0.17,-0.19), (0.13,-0.39), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (0.11,-0.03), (-0.17,0.46), (-0.17,-0.19), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (0.19,-0.17), (0.32,0.09), (0.23,-0.24), (0.18,0.01), (2.0,3.0), (2.0,3.0), (2.0,3.0), (2.0,3.0)/;
      DATA              ((CTRUE5(I,J,2),I=1,8),J=1,5)/(0.1,0.1), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (-0.16,-0.37), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (-0.17,-0.19), (8.0,9.0), (0.13,-0.39), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (0.11,-0.03), (3.0,6.0), (-0.17,0.46), (4.0,7.0), (-0.17,-0.19), (7.0,2.0), (7.0,2.0), (7.0,2.0), (0.19,-0.17), (5.0,8.0), (0.32,0.09), (6.0,9.0), (0.23,-0.24), (8.0,3.0), (0.18,0.01), (9.0,4.0)/;
      DATA              ((CTRUE6(I,J,1),I=1,8),J=1,5)/(0.1,0.1), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (1.0,2.0), (0.09,-0.12), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (3.0,4.0), (0.03,-0.09), (0.15,-0.03), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (5.0,6.0), (0.03,0.03), (-0.18,0.03), (0.03,-0.09), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (7.0,8.0), (0.09,0.03), (0.03,0.12), (0.12,0.03), (0.03,0.06), (2.0,3.0), (2.0,3.0), (2.0,3.0), (2.0,3.0)/;
      DATA              ((CTRUE6(I,J,2),I=1,8),J=1,5)/(0.1,0.1), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (4.0,5.0), (0.09,-0.12), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (6.0,7.0), (0.03,-0.09), (8.0,9.0), (0.15,-0.03), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (2.0,5.0), (0.03,0.03), (3.0,6.0), (-0.18,0.03), (4.0,7.0), (0.03,-0.09), (7.0,2.0), (7.0,2.0), (7.0,2.0), (0.09,0.03), (5.0,8.0), (0.03,0.12), (6.0,9.0), (0.12,0.03), (8.0,3.0), (0.03,0.06), (9.0,4.0)/;
      const ITRUE3 = [0, 1, 2, 2, 2];
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 60
         for (NP1 = 1; NP1 <= 5; NP1++) { // 40
            N = NP1 - 1;
            LEN = 2*max(N,1);
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               CX[I] = CV(I,NP1,INCX);
            } // 20
            if (ICASE == 6) {
               // .. SCNRM2TEST ..
               stest1(SCNRM2TEST(N,CX,INCX),STRUE2(NP1), STRUE2(NP1), SFAC);
            } else if (ICASE == 7) {
               // .. SCASUMTEST ..
               stest1(SCASUMTEST(N,CX,INCX),STRUE4(NP1), STRUE4(NP1),SFAC);
            } else if (ICASE == 8) {
               // .. CSCAL ..
               cscal(N,CA,CX,INCX);
               ctest(LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX), SFAC);
            } else if (ICASE == 9) {
               // .. CSSCALTEST ..
               csscaltest(N,SA,CX,INCX);
               ctest(LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX), SFAC);
            } else if (ICASE == 10) {
               // .. ICAMAXTEST ..
               itest1(ICAMAXTEST(N,CX,INCX),ITRUE3(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1';
               STOP;
            }

         } // 40
      } // 60

      INCX = 1;
      if (ICASE == 8) {
         // CSCAL
         // Add a test for alpha equal to zero.
         CA = (0.0,0.0);
         for (I = 1; I <= 5; I++) { // 80
            MWPCT[I] = (0.0,0.0);
            MWPCS[I] = (1.0,1.0);
         } // 80
         cscal(5,CA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      } else if (ICASE == 9) {
         // CSSCALTEST
         // Add a test for alpha equal to zero.
         SA = 0.0;
         for (I = 1; I <= 5; I++) { // 100
            MWPCT[I] = (0.0,0.0);
            MWPCS[I] = (1.0,1.0);
         } // 100
         csscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to one.
         SA = 1.0;
         for (I = 1; I <= 5; I++) { // 120
            MWPCT[I] = CX(I);
            MWPCS[I] = CX(I);
         } // 120
         csscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to minus one.
         SA = -1.0;
         for (I = 1; I <= 5; I++) { // 140
            MWPCT[I] = -CX(I);
            MWPCS[I] = -CX(I);
         } // 140
         csscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      }
      return;
      }
      void check2(SFAC) {
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double              SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      Complex           CA,CTEMP;
      int               I, J, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      Complex           CDOT(1), CSIZE1(4), CSIZE2(7,2), CSIZE3(14), CT10X(7,4,4), CT10Y(7,4,4), CT6(4,4), CT7(4,4), CT8(7,4,4), CX(7), CX1(7), CY(7), CY1(7);
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      // EXTERNAL CDOTCTEST, CDOTUTEST
      // .. External Subroutines ..
      // EXTERNAL CAXPYTEST, CCOPYTEST, CSWAPTEST, CTEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      const CA = [(0.4,-0.7)];
      const INCXS = [1, 2, -2, -1];
      const INCYS = [1, -2, 1, -2];
      const LENS = [1, 1, 2, 4, 1, 1, 3, 7];
      const NS = [0, 1, 2, 4];
      const CX1 = [(0.7,-0.8), (-0.4,-0.7), (-0.1,-0.9), (0.2,-0.8), (-0.9,-0.4), (0.1,0.4), (-0.6,0.6)];
      const CY1 = [(0.6,-0.6), (-0.9,0.5), (0.7,-0.6), (0.1,-0.5), (-0.1,-0.2), (-0.5,-0.3), (0.8,-0.7)];
      DATA              ((CT8(I,J,1),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (-1.55,0.5), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (-1.55,0.5), (0.03,-0.89), (-0.38,-0.96), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT8(I,J,2),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.07,-0.89), (-0.9,0.5), (0.42,-1.41), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.78,0.06), (-0.9,0.5), (0.06,-0.13), (0.1,-0.5), (-0.77,-0.49), (-0.5,-0.3), (0.52,-1.51)/;
      DATA              ((CT8(I,J,3),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.07,-0.89), (-1.18,-0.31), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.78,0.06), (-1.54,0.97), (0.03,-0.89), (-0.18,-1.31), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT8(I,J,4),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (-0.9,0.5), (0.05,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.32,-1.41), (-0.9,0.5), (0.05,-0.6), (0.1,-0.5), (-0.77,-0.49), (-0.5,-0.3), (0.32,-1.16)/;
      const CT7 = [(0.0,0.0), (-0.06,-0.90), (0.65,-0.47), (-0.34,-1.22), (0.0,0.0), (-0.06,-0.90), (-0.59,-1.46), (-1.04,-0.04), (0.0,0.0), (-0.06,-0.90), (-0.83,0.59), (0.07,-0.37), (0.0,0.0), (-0.06,-0.90), (-0.76,-1.15), (-1.33,-1.82)];
      const CT6 = [(0.0,0.0), (0.90,0.06), (0.91,-0.77), (1.80,-0.10), (0.0,0.0), (0.90,0.06), (1.45,0.74), (0.20,0.90), (0.0,0.0), (0.90,0.06), (-0.55,0.23), (0.83,-0.39), (0.0,0.0), (0.90,0.06), (1.04,0.79), (1.95,1.22)];
      DATA              ((CT10X(I,J,1),I=1,7),J=1,4)/(0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (-0.9,0.5), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (-0.9,0.5), (0.7,-0.6), (0.1,-0.5), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT10X(I,J,2),I=1,7),J=1,4)/(0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.6), (-0.4,-0.7), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.8,-0.7), (-0.4,-0.7), (-0.1,-0.2), (0.2,-0.8), (0.7,-0.6), (0.1,0.4), (0.6,-0.6)/;
      DATA              ((CT10X(I,J,3),I=1,7),J=1,4)/(0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.9,0.5), (-0.4,-0.7), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.1,-0.5), (-0.4,-0.7), (0.7,-0.6), (0.2,-0.8), (-0.9,0.5), (0.1,0.4), (0.6,-0.6)/;
      DATA              ((CT10X(I,J,4),I=1,7),J=1,4)/(0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.7,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.6,-0.6), (0.7,-0.6), (-0.1,-0.2), (0.8,-0.7), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT10Y(I,J,1),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (-0.4,-0.7), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (-0.4,-0.7), (-0.1,-0.9), (0.2,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT10Y(I,J,2),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.1,-0.9), (-0.9,0.5), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.6,0.6), (-0.9,0.5), (-0.9,-0.4), (0.1,-0.5), (-0.1,-0.9), (-0.5,-0.3), (0.7,-0.8)/;
      DATA              ((CT10Y(I,J,3),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.1,-0.9), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (-0.6,0.6), (-0.9,-0.4), (-0.1,-0.9), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0)/;
      DATA              ((CT10Y(I,J,4),I=1,7),J=1,4)/(0.6,-0.6), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (-0.9,0.5), (-0.4,-0.7), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.7,-0.8), (-0.9,0.5), (-0.4,-0.7), (0.1,-0.5), (-0.1,-0.9), (-0.5,-0.3), (0.2,-0.8)/;
      const CSIZE1 = [(0.0,0.0), (0.9,0.9), (1.63,1.73), (2.90,2.78)];
      const CSIZE3 = [(0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (1.17,1.17), (1.17,1.17), (1.17,1.17), (1.17,1.17), (1.17,1.17), (1.17,1.17), (1.17,1.17)];
      const CSIZE2 = [(0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (0.0,0.0), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54)];
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
            // .. initialize all argument arrays ..
            for (I = 1; I <= 7; I++) { // 20
               CX[I] = CX1(I);
               CY[I] = CY1(I);
            } // 20
            if (ICASE == 1) {
               // .. CDOTCTEST ..
               cdotctest(N,CX,INCX,CY,INCY,CTEMP);
               CDOT[1] = CTEMP;
               ctest(1,CDOT,CT6(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE == 2) {
               // .. CDOTUTEST ..
               cdotutest(N,CX,INCX,CY,INCY,CTEMP);
               CDOT[1] = CTEMP;
               ctest(1,CDOT,CT7(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE == 3) {
               // .. CAXPYTEST ..
               caxpytest(N,CA,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC);
            } else if (ICASE == 4) {
               // .. CCOPYTEST ..
               ccopytest(N,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0);
            } else if (ICASE == 5) {
               // .. CSWAPTEST ..
               cswaptest(N,CX,INCX,CY,INCY);
               ctest(LENX,CX,CT10X(1,KN,KI),CSIZE3,1.0);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2';
               STOP;
            }

         } // 40
      } // 60
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
      const            NOUT=6;
      // .. Scalar Arguments ..
      double             SFAC;
      int              LEN;
      // .. Array Arguments ..
      double             SCOMP(LEN), SSIZE(LEN), STRUE(LEN);
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      double             SD;
      int              I;
      // .. External Functions ..
      //- REAL             SDIFF;
      // EXTERNAL SDIFF
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I);
         if (SDIFF((SSIZE(I)).abs()+(SFAC*SD).abs(),(SSIZE(I))).abs() == 0.0) GO TO 40;

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
      void stest1(SCOMP1,STRUE1,SSIZE,SFAC) {
      // ************************* STEST1 *****************************

      // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
      // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
      // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      double              SCOMP1, SFAC, STRUE1;
      // .. Array Arguments ..
      double              SSIZE(*);
      // .. Local Arrays ..
      double              SCOMP(1), STRUE(1);
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Executable Statements ..

      SCOMP[1] = SCOMP1;
      STRUE[1] = STRUE1;
      stest(1,SCOMP,STRUE,SSIZE,SFAC);

      return;
      }
      double sdiff(SA,SB) {
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      double                            SA, SB;
      // .. Executable Statements ..
      SDIFF = SA - SB;
      return;
      }
      void ctest(LEN,CCOMP,CTRUE,CSIZE,SFAC) {
      // **************************** CTEST *****************************

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      double             SFAC;
      int              LEN;
      // .. Array Arguments ..
      Complex          CCOMP(LEN), CSIZE(LEN), CTRUE(LEN);
      // .. Local Scalars ..
      int              I;
      // .. Local Arrays ..
      double             SCOMP(20), SSIZE(20), STRUE(20);
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
      // .. Executable Statements ..
      for (I = 1; I <= LEN; I++) { // 20
         SCOMP[2*I-1] = double(CCOMP(I));
         SCOMP[2*I] = AIMAG(CCOMP(I));
         STRUE[2*I-1] = double(CTRUE(I));
         STRUE[2*I] = AIMAG(CTRUE(I));
         SSIZE[2*I-1] = double(CSIZE(I));
         SSIZE[2*I] = AIMAG(CSIZE(I));
      } // 20

      stest(2*LEN,SCOMP,STRUE,SSIZE,SFAC);
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
