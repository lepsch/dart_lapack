void main() {      // Test program for the COMPLEX*16 Level 1 CBLAS.
      // Based upon the original CBLAS test routine together with:
      // F06GAF Example Program Text
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      double           SFAC;
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK1, CHECK2, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA             SFAC/9.765625D-4/
      // .. Executable Statements ..
      WRITE (NOUT,99999)
      for (IC = 1; IC <= 10; IC++) { // 20
         ICASE = IC
         header();

         // Initialize PASS, INCX, INCY, and MODE for a new case.
         // The value 9999 for INCX, INCY or MODE will appear in the
         // detailed  output, if any, for cases that do not involve
         // these parameters.

         PASS = .TRUE.
         INCX = 9999
         INCY = 9999
         MODE = 9999
         if (ICASE.LE.5) {
            check2(SFAC);
         } else if (ICASE.GE.6) {
            check1(SFAC);
         }
         // -- Print
         if (PASS) WRITE (NOUT,99998);
      } // 20
      STOP

99999 FORMAT (' Complex CBLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')
      }
      SUBROUTINE HEADER
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
      DATA             L(1)/'CBLAS_ZDOTC'/
      DATA             L(2)/'CBLAS_ZDOTU'/
      DATA             L(3)/'CBLAS_ZAXPY'/
      DATA             L(4)/'CBLAS_ZCOPY'/
      DATA             L(5)/'CBLAS_ZSWAP'/
      DATA             L(6)/'CBLAS_DZNRM2'/
      DATA             L(7)/'CBLAS_DZASUM'/
      DATA             L(8)/'CBLAS_ZSCAL'/
      DATA             L(9)/'CBLAS_ZDSCAL'/
      DATA             L(10)/'CBLAS_IZAMAX'/
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE)
      RETURN

99999 FORMAT (/' Test of subprogram number',I3,9X,A15)
      }
      SUBROUTINE CHECK1(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      COMPLEX*16        CA
      double            SA;
      int               I, J, LEN, NP1;
      // .. Local Arrays ..
      COMPLEX*16        CTRUE5(8,5,2), CTRUE6(8,5,2), CV(8,5,2), CX(8), MWPCS(5), MWPCT(5)
      double            STRUE2(5), STRUE4(5);
      int               ITRUE3(5);
      // .. External Functions ..
      double            DZASUMTEST, DZNRM2TEST;
      int               IZAMAXTEST;
      // EXTERNAL DZASUMTEST, DZNRM2TEST, IZAMAXTEST
      // .. External Subroutines ..
      // EXTERNAL ZSCALTEST, ZDSCALTEST, CTEST, ITEST1, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              SA, CA/0.3D0, (0.4D0,-0.7D0)/
      DATA              ((CV(I,J,1),I=1,8),J=1,5)/(0.1D0,0.1D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (0.3D0,-0.4D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (0.1D0,-0.3D0), (0.5D0,-0.1D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (0.1D0,0.1D0), (-0.6D0,0.1D0), (0.1D0,-0.3D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (0.3D0,0.1D0), (0.1D0,0.4D0), (0.4D0,0.1D0), (0.1D0,0.2D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0)/
      DATA              ((CV(I,J,2),I=1,8),J=1,5)/(0.1D0,0.1D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (0.3D0,-0.4D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (0.1D0,-0.3D0), (8.0D0,9.0D0), (0.5D0,-0.1D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (0.1D0,0.1D0), (3.0D0,6.0D0), (-0.6D0,0.1D0), (4.0D0,7.0D0), (0.1D0,-0.3D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (0.3D0,0.1D0), (5.0D0,8.0D0), (0.1D0,0.4D0), (6.0D0,9.0D0), (0.4D0,0.1D0), (8.0D0,3.0D0), (0.1D0,0.2D0), (9.0D0,4.0D0)/
      DATA              STRUE2/0.0D0, 0.5D0, 0.6D0, 0.7D0, 0.7D0/
      DATA              STRUE4/0.0D0, 0.7D0, 1.0D0, 1.3D0, 1.7D0/
      DATA              ((CTRUE5(I,J,1),I=1,8),J=1,5)/(0.1D0,0.1D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (-0.16D0,-0.37D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (-0.17D0,-0.19D0), (0.13D0,-0.39D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (0.11D0,-0.03D0), (-0.17D0,0.46D0), (-0.17D0,-0.19D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (0.19D0,-0.17D0), (0.32D0,0.09D0), (0.23D0,-0.24D0), (0.18D0,0.01D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0)/
      DATA              ((CTRUE5(I,J,2),I=1,8),J=1,5)/(0.1D0,0.1D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (-0.16D0,-0.37D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (-0.17D0,-0.19D0), (8.0D0,9.0D0), (0.13D0,-0.39D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (0.11D0,-0.03D0), (3.0D0,6.0D0), (-0.17D0,0.46D0), (4.0D0,7.0D0), (-0.17D0,-0.19D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (0.19D0,-0.17D0), (5.0D0,8.0D0), (0.32D0,0.09D0), (6.0D0,9.0D0), (0.23D0,-0.24D0), (8.0D0,3.0D0), (0.18D0,0.01D0), (9.0D0,4.0D0)/
      DATA              ((CTRUE6(I,J,1),I=1,8),J=1,5)/(0.1D0,0.1D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (1.0D0,2.0D0), (0.09D0,-0.12D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (3.0D0,4.0D0), (0.03D0,-0.09D0), (0.15D0,-0.03D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (5.0D0,6.0D0), (0.03D0,0.03D0), (-0.18D0,0.03D0), (0.03D0,-0.09D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (7.0D0,8.0D0), (0.09D0,0.03D0), (0.03D0,0.12D0), (0.12D0,0.03D0), (0.03D0,0.06D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0), (2.0D0,3.0D0)/
      DATA              ((CTRUE6(I,J,2),I=1,8),J=1,5)/(0.1D0,0.1D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (4.0D0,5.0D0), (0.09D0,-0.12D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (6.0D0,7.0D0), (0.03D0,-0.09D0), (8.0D0,9.0D0), (0.15D0,-0.03D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (2.0D0,5.0D0), (0.03D0,0.03D0), (3.0D0,6.0D0), (-0.18D0,0.03D0), (4.0D0,7.0D0), (0.03D0,-0.09D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (7.0D0,2.0D0), (0.09D0,0.03D0), (5.0D0,8.0D0), (0.03D0,0.12D0), (6.0D0,9.0D0), (0.12D0,0.03D0), (8.0D0,3.0D0), (0.03D0,0.06D0), (9.0D0,4.0D0)/
      DATA              ITRUE3/0, 1, 2, 2, 2/
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 60
         for (NP1 = 1; NP1 <= 5; NP1++) { // 40
            N = NP1 - 1
            LEN = 2*MAX(N,1)
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               CX(I) = CV(I,NP1,INCX)
            } // 20
            if (ICASE.EQ.6) {
               // .. DZNRM2TEST ..
               stest1(DZNRM2TEST(N,CX,INCX),STRUE2(NP1), STRUE2(NP1),SFAC);
            } else if (ICASE.EQ.7) {
               // .. DZASUMTEST ..
               stest1(DZASUMTEST(N,CX,INCX),STRUE4(NP1), STRUE4(NP1),SFAC);
            } else if (ICASE.EQ.8) {
               // .. ZSCALTEST ..
               zscaltest(N,CA,CX,INCX);
               ctest(LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX), SFAC);
            } else if (ICASE.EQ.9) {
               // .. ZDSCALTEST ..
               zdscaltest(N,SA,CX,INCX);
               ctest(LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX), SFAC);
            } else if (ICASE.EQ.10) {
               // .. IZAMAXTEST ..
               itest1(IZAMAXTEST(N,CX,INCX),ITRUE3(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            }

         } // 40
      } // 60

      INCX = 1
      if (ICASE.EQ.8) {
         // ZSCALTEST
         // Add a test for alpha equal to zero.
         CA = (0.0D0,0.0D0)
         for (I = 1; I <= 5; I++) { // 80
            MWPCT(I) = (0.0D0,0.0D0)
            MWPCS(I) = (1.0D0,1.0D0)
         } // 80
         zscaltest(5,CA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      } else if (ICASE.EQ.9) {
         // ZDSCALTEST
         // Add a test for alpha equal to zero.
         SA = 0.0D0
         for (I = 1; I <= 5; I++) { // 100
            MWPCT(I) = (0.0D0,0.0D0)
            MWPCS(I) = (1.0D0,1.0D0)
         } // 100
         zdscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to one.
         SA = 1.0D0
         for (I = 1; I <= 5; I++) { // 120
            MWPCT(I) = CX(I)
            MWPCS(I) = CX(I)
         } // 120
         zdscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to minus one.
         SA = -1.0D0
         for (I = 1; I <= 5; I++) { // 140
            MWPCT(I) = -CX(I)
            MWPCS(I) = -CX(I)
         } // 140
         zdscaltest(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      }
      RETURN
      }
      SUBROUTINE CHECK2(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      double            SFAC;
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      COMPLEX*16        CA,ZTEMP
      int               I, J, KI, KN, KSIZE, LENX, LENY, MX, MY;
      // .. Local Arrays ..
      COMPLEX*16        CDOT(1), CSIZE1(4), CSIZE2(7,2), CSIZE3(14), CT10X(7,4,4), CT10Y(7,4,4), CT6(4,4), CT7(4,4), CT8(7,4,4), CX(7), CX1(7), CY(7), CY1(7)
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      // EXTERNAL ZDOTCTEST, ZDOTUTEST
      // .. External Subroutines ..
      // EXTERNAL ZAXPYTEST, ZCOPYTEST, ZSWAPTEST, CTEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              CA/(0.4D0,-0.7D0)/
      DATA              INCXS/1, 2, -2, -1/
      DATA              INCYS/1, -2, 1, -2/
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/
      DATA              NS/0, 1, 2, 4/
      DATA              CX1/(0.7D0,-0.8D0), (-0.4D0,-0.7D0), (-0.1D0,-0.9D0), (0.2D0,-0.8D0), (-0.9D0,-0.4D0), (0.1D0,0.4D0), (-0.6D0,0.6D0)/
      DATA              CY1/(0.6D0,-0.6D0), (-0.9D0,0.5D0), (0.7D0,-0.6D0), (0.1D0,-0.5D0), (-0.1D0,-0.2D0), (-0.5D0,-0.3D0), (0.8D0,-0.7D0)/
      DATA              ((CT8(I,J,1),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (-1.55D0,0.5D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (-1.55D0,0.5D0), (0.03D0,-0.89D0), (-0.38D0,-0.96D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT8(I,J,2),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.07D0,-0.89D0), (-0.9D0,0.5D0), (0.42D0,-1.41D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.78D0,0.06D0), (-0.9D0,0.5D0), (0.06D0,-0.13D0), (0.1D0,-0.5D0), (-0.77D0,-0.49D0), (-0.5D0,-0.3D0), (0.52D0,-1.51D0)/
      DATA              ((CT8(I,J,3),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.07D0,-0.89D0), (-1.18D0,-0.31D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.78D0,0.06D0), (-1.54D0,0.97D0), (0.03D0,-0.89D0), (-0.18D0,-1.31D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT8(I,J,4),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (-0.9D0,0.5D0), (0.05D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.32D0,-1.41D0), (-0.9D0,0.5D0), (0.05D0,-0.6D0), (0.1D0,-0.5D0), (-0.77D0,-0.49D0), (-0.5D0,-0.3D0), (0.32D0,-1.16D0)/
      DATA              CT7/(0.0D0,0.0D0), (-0.06D0,-0.90D0), (0.65D0,-0.47D0), (-0.34D0,-1.22D0), (0.0D0,0.0D0), (-0.06D0,-0.90D0), (-0.59D0,-1.46D0), (-1.04D0,-0.04D0), (0.0D0,0.0D0), (-0.06D0,-0.90D0), (-0.83D0,0.59D0), (0.07D0,-0.37D0), (0.0D0,0.0D0), (-0.06D0,-0.90D0), (-0.76D0,-1.15D0), (-1.33D0,-1.82D0)/
      DATA              CT6/(0.0D0,0.0D0), (0.90D0,0.06D0), (0.91D0,-0.77D0), (1.80D0,-0.10D0), (0.0D0,0.0D0), (0.90D0,0.06D0), (1.45D0,0.74D0), (0.20D0,0.90D0), (0.0D0,0.0D0), (0.90D0,0.06D0), (-0.55D0,0.23D0), (0.83D0,-0.39D0), (0.0D0,0.0D0), (0.90D0,0.06D0), (1.04D0,0.79D0), (1.95D0,1.22D0)/
      DATA              ((CT10X(I,J,1),I=1,7),J=1,4)/(0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (-0.9D0,0.5D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (-0.9D0,0.5D0), (0.7D0,-0.6D0), (0.1D0,-0.5D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT10X(I,J,2),I=1,7),J=1,4)/(0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.6D0), (-0.4D0,-0.7D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.8D0,-0.7D0), (-0.4D0,-0.7D0), (-0.1D0,-0.2D0), (0.2D0,-0.8D0), (0.7D0,-0.6D0), (0.1D0,0.4D0), (0.6D0,-0.6D0)/
      DATA              ((CT10X(I,J,3),I=1,7),J=1,4)/(0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.9D0,0.5D0), (-0.4D0,-0.7D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.1D0,-0.5D0), (-0.4D0,-0.7D0), (0.7D0,-0.6D0), (0.2D0,-0.8D0), (-0.9D0,0.5D0), (0.1D0,0.4D0), (0.6D0,-0.6D0)/
      DATA              ((CT10X(I,J,4),I=1,7),J=1,4)/(0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.7D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.6D0,-0.6D0), (0.7D0,-0.6D0), (-0.1D0,-0.2D0), (0.8D0,-0.7D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT10Y(I,J,1),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (-0.4D0,-0.7D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (-0.4D0,-0.7D0), (-0.1D0,-0.9D0), (0.2D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT10Y(I,J,2),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.1D0,-0.9D0), (-0.9D0,0.5D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.6D0,0.6D0), (-0.9D0,0.5D0), (-0.9D0,-0.4D0), (0.1D0,-0.5D0), (-0.1D0,-0.9D0), (-0.5D0,-0.3D0), (0.7D0,-0.8D0)/
      DATA              ((CT10Y(I,J,3),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.1D0,-0.9D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (-0.6D0,0.6D0), (-0.9D0,-0.4D0), (-0.1D0,-0.9D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0)/
      DATA              ((CT10Y(I,J,4),I=1,7),J=1,4)/(0.6D0,-0.6D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (-0.9D0,0.5D0), (-0.4D0,-0.7D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.7D0,-0.8D0), (-0.9D0,0.5D0), (-0.4D0,-0.7D0), (0.1D0,-0.5D0), (-0.1D0,-0.9D0), (-0.5D0,-0.3D0), (0.2D0,-0.8D0)/
      DATA              CSIZE1/(0.0D0,0.0D0), (0.9D0,0.9D0), (1.63D0,1.73D0), (2.90D0,2.78D0)/
      DATA              CSIZE3/(0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (1.17D0,1.17D0), (1.17D0,1.17D0), (1.17D0,1.17D0), (1.17D0,1.17D0), (1.17D0,1.17D0), (1.17D0,1.17D0), (1.17D0,1.17D0)/
      DATA              CSIZE2/(0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (0.0D0,0.0D0), (1.54D0,1.54D0), (1.54D0,1.54D0), (1.54D0,1.54D0), (1.54D0,1.54D0), (1.54D0,1.54D0), (1.54D0,1.54D0), (1.54D0,1.54D0)/
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
            // .. initialize all argument arrays ..
            for (I = 1; I <= 7; I++) { // 20
               CX(I) = CX1(I)
               CY(I) = CY1(I)
            } // 20
            if (ICASE.EQ.1) {
               // .. ZDOTCTEST ..
               zdotctest(N,CX,INCX,CY,INCY,ZTEMP);
               CDOT(1) = ZTEMP
               ctest(1,CDOT,CT6(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE.EQ.2) {
               // .. ZDOTUTEST ..
               zdotutest(N,CX,INCX,CY,INCY,ZTEMP);
               CDOT(1) = ZTEMP
               ctest(1,CDOT,CT7(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE.EQ.3) {
               // .. ZAXPYTEST ..
               zaxpytest(N,CA,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC);
            } else if (ICASE.EQ.4) {
               // .. ZCOPYTEST ..
               zcopytest(N,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0D0);
            } else if (ICASE.EQ.5) {
               // .. ZSWAPTEST ..
               zswaptest(N,CX,INCX,CY,INCY);
               ctest(LENX,CX,CT10X(1,KN,KI),CSIZE3,1.0D0);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0D0);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            }

         } // 40
      } // 60
      RETURN
      }
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC)
      // ********************************* STEST **************************

      // THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
      // SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
      // NEGLIGIBLE.

      // C. L. LAWSON, JPL, 1974 DEC 10

      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalar Arguments ..
      double           SFAC;
      int              LEN;
      // .. Array Arguments ..
      double           SCOMP(LEN), SSIZE(LEN), STRUE(LEN);
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
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
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I)
         IF (SDIFF(ABS(SSIZE(I))+ABS(SFAC*SD),ABS(SSIZE(I))).EQ.0.0D0) GO TO 40

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = .FALSE.
         WRITE (NOUT,99999)
         WRITE (NOUT,99998)
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, I, SCOMP(I), STRUE(I), SD, SSIZE(I)
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE  I                            ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,3I5,I3,2D36.8,2D12.4)
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
      }
      double           FUNCTION SDIFF(SA,SB);
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      double                          SA, SB;
      // .. Executable Statements ..
      SDIFF = SA - SB
      RETURN
      }
      SUBROUTINE CTEST(LEN,CCOMP,CTRUE,CSIZE,SFAC)
      // **************************** CTEST *****************************

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      double           SFAC;
      int              LEN;
      // .. Array Arguments ..
      COMPLEX*16       CCOMP(LEN), CSIZE(LEN), CTRUE(LEN)
      // .. Local Scalars ..
      int              I;
      // .. Local Arrays ..
      double           SCOMP(20), SSIZE(20), STRUE(20);
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Intrinsic Functions ..
      // INTRINSIC DIMAG, DBLE
      // .. Executable Statements ..
      for (I = 1; I <= LEN; I++) { // 20
         SCOMP(2*I-1) = DBLE(CCOMP(I))
         SCOMP(2*I) = DIMAG(CCOMP(I))
         STRUE(2*I-1) = DBLE(CTRUE(I))
         STRUE(2*I) = DIMAG(CTRUE(I))
         SSIZE(2*I-1) = DBLE(CSIZE(I))
         SSIZE(2*I) = DIMAG(CSIZE(I))
      } // 20

      stest(2*LEN,SCOMP,STRUE,SSIZE,SFAC);
      RETURN
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
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      int               ID;
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..
      if (ICOMP.EQ.ITRUE) GO TO 40;

                             // HERE ICOMP IS NOT EQUAL TO ITRUE.

      if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
      PASS = .FALSE.
      WRITE (NOUT,99999)
      WRITE (NOUT,99998)
   20 ID = ICOMP - ITRUE
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, ICOMP, ITRUE, ID
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE                               ', ' COMP                                TRUE     DIFFERENCE', /1X)
99997 FORMAT (1X,I4,I3,3I5,2I36,I12)
      }
