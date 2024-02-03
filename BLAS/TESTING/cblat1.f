void main() {
*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SFAC
      int              IC;
      // .. External Subroutines ..
      // EXTERNAL CHECK1, CHECK2, HEADER
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA             SFAC/9.765625E-4/
      // .. Executable Statements ..
      WRITE (NOUT,99999)
      for (IC = 1; IC <= 10; IC++) { // 20
         ICASE = IC
         header();

         // Initialize PASS, INCX, INCY, and MODE for a new case.
         // The value 9999 for INCX, INCY or MODE will appear in the
         // detailed  output, if any, for cases that do not involve
         // these parameters.

         PASS = true;
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

99999 FORMAT (' Complex BLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')

      // End of CBLAT1

      }
      SUBROUTINE HEADER
      // .. Parameters ..
      int              NOUT;
      const            NOUT=6;
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Arrays ..
      String           L(10);
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA             L(1)/'CDOTC '/
      DATA             L(2)/'CDOTU '/
      DATA             L(3)/'CAXPY '/
      DATA             L(4)/'CCOPY '/
      DATA             L(5)/'CSWAP '/
      DATA             L(6)/'SCNRM2'/
      DATA             L(7)/'SCASUM'/
      DATA             L(8)/'CSCAL '/
      DATA             L(9)/'CSSCAL'/
      DATA             L(10)/'ICAMAX'/
      // .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE)
      RETURN

99999 FORMAT (/' Test of subprogram number',I3,12X,A6)

      // End of HEADER

      }
      SUBROUTINE CHECK1(SFAC)
      // .. Parameters ..
      int               NOUT;
      REAL              THRESH
      const             NOUT=6, THRESH=10.0E0;
      // .. Scalar Arguments ..
      REAL              SFAC
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      COMPLEX           CA
      REAL              SA
      int               I, IX, J, LEN, NP1;
      // .. Local Arrays ..
      COMPLEX           CTRUE5(8,5,2), CTRUE6(8,5,2), CV(8,5,2), CVR(8), CX(8), CXR(15), MWPCS(5), MWPCT(5)
      REAL              STRUE2(5), STRUE4(5)
      int               ITRUE3(5), ITRUEC(5);
      // .. External Functions ..
      REAL              SCASUM, SCNRM2
      int               ICAMAX;
      // EXTERNAL SCASUM, SCNRM2, ICAMAX
      // .. External Subroutines ..
      // EXTERNAL CB1NRM2, CSCAL, CSSCAL, CTEST, ITEST1, STEST1
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              SA, CA/0.3E0, (0.4E0,-0.7E0)/
      DATA              ((CV(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (0.3E0,-0.4E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (0.1E0,-0.3E0), (0.5E0,-0.1E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (0.1E0,0.1E0), (-0.6E0,0.1E0), (0.1E0,-0.3E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (0.3E0,0.1E0), (0.5E0,0.0E0), (0.0E0,0.5E0), (0.0E0,0.2E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CV(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (0.3E0,-0.4E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (0.1E0,-0.3E0), (8.0E0,9.0E0), (0.5E0,-0.1E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (0.1E0,0.1E0), (3.0E0,6.0E0), (-0.6E0,0.1E0), (4.0E0,7.0E0), (0.1E0,-0.3E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (0.3E0,0.1E0), (5.0E0,8.0E0), (0.5E0,0.0E0), (6.0E0,9.0E0), (0.0E0,0.5E0), (8.0E0,3.0E0), (0.0E0,0.2E0), (9.0E0,4.0E0)/
      DATA              CVR/(8.0E0,8.0E0), (-7.0E0,-7.0E0), (9.0E0,9.0E0), (5.0E0,5.0E0), (9.0E0,9.0E0), (8.0E0,8.0E0), (7.0E0,7.0E0), (7.0E0,7.0E0)/
      DATA              STRUE2/0.0E0, 0.5E0, 0.6E0, 0.7E0, 0.8E0/
      DATA              STRUE4/0.0E0, 0.7E0, 1.0E0, 1.3E0, 1.6E0/
      DATA              ((CTRUE5(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (-0.16E0,-0.37E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (-0.17E0,-0.19E0), (0.13E0,-0.39E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (0.11E0,-0.03E0), (-0.17E0,0.46E0), (-0.17E0,-0.19E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (0.19E0,-0.17E0), (0.20E0,-0.35E0), (0.35E0,0.20E0), (0.14E0,0.08E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CTRUE5(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (-0.16E0,-0.37E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (-0.17E0,-0.19E0), (8.0E0,9.0E0), (0.13E0,-0.39E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (0.11E0,-0.03E0), (3.0E0,6.0E0), (-0.17E0,0.46E0), (4.0E0,7.0E0), (-0.17E0,-0.19E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (0.19E0,-0.17E0), (5.0E0,8.0E0), (0.20E0,-0.35E0), (6.0E0,9.0E0), (0.35E0,0.20E0), (8.0E0,3.0E0), (0.14E0,0.08E0), (9.0E0,4.0E0)/
      DATA              ((CTRUE6(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0), (0.09E0,-0.12E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0), (0.03E0,-0.09E0), (0.15E0,-0.03E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0), (0.03E0,0.03E0), (-0.18E0,0.03E0), (0.03E0,-0.09E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0), (0.09E0,0.03E0), (0.15E0,0.00E0), (0.00E0,0.15E0), (0.00E0,0.06E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CTRUE6(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0), (0.09E0,-0.12E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0), (0.03E0,-0.09E0), (8.0E0,9.0E0), (0.15E0,-0.03E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0), (0.03E0,0.03E0), (3.0E0,6.0E0), (-0.18E0,0.03E0), (4.0E0,7.0E0), (0.03E0,-0.09E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (7.0E0,2.0E0), (0.09E0,0.03E0), (5.0E0,8.0E0), (0.15E0,0.00E0), (6.0E0,9.0E0), (0.00E0,0.15E0), (8.0E0,3.0E0), (0.00E0,0.06E0), (9.0E0,4.0E0)/
      DATA              ITRUE3/0, 1, 2, 2, 2/
      DATA              ITRUEC/0, 1, 1, 1, 1/
      // .. Executable Statements ..
      for (INCX = 1; INCX <= 2; INCX++) { // 60
         for (NP1 = 1; NP1 <= 5; NP1++) { // 40
            N = NP1 - 1
            LEN = 2*MAX(N,1)
            // .. Set vector arguments ..
            for (I = 1; I <= LEN; I++) { // 20
               CX(I) = CV(I,NP1,INCX)
            } // 20
            if (ICASE == 6) {
               // .. SCNRM2 ..
               // Test scaling when some entries are tiny or huge
               cb1nrm2(N,(INCX-2)*2,THRESH);
               cb1nrm2(N,INCX,THRESH);
               // Test with hardcoded mid range entries
               stest1(SCNRM2(N,CX,INCX),STRUE2(NP1),STRUE2(NP1), SFAC);
            } else if (ICASE == 7) {
               // .. SCASUM ..
               stest1(SCASUM(N,CX,INCX),STRUE4(NP1),STRUE4(NP1), SFAC);
            } else if (ICASE == 8) {
               // .. CSCAL ..
               cscal(N,CA,CX,INCX);
               ctest(LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX), SFAC);
            } else if (ICASE == 9) {
               // .. CSSCAL ..
               csscal(N,SA,CX,INCX);
               ctest(LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX), SFAC);
            } else if (ICASE == 10) {
               // .. ICAMAX ..
               itest1(ICAMAX(N,CX,INCX),ITRUE3(NP1));
               for (I = 1; I <= LEN; I++) { // 160
                  CX(I) = (42.0E0,43.0E0)
               } // 160
               itest1(ICAMAX(N,CX,INCX),ITRUEC(NP1));
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            }

         } // 40
         if (ICASE == 10) {
            N = 8
            IX = 1
            for (I = 1; I <= N; I++) { // 180
               CXR(IX) = CVR(I)
               IX = IX + INCX
            } // 180
            itest1(ICAMAX(N,CXR,INCX),3);
         }
      } // 60

      INCX = 1
      if (ICASE == 8) {
         // CSCAL
         // Add a test for alpha equal to zero.
         CA = (0.0E0,0.0E0)
         for (I = 1; I <= 5; I++) { // 80
            MWPCT(I) = (0.0E0,0.0E0)
            MWPCS(I) = (1.0E0,1.0E0)
         } // 80
         cscal(5,CA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      } else if (ICASE == 9) {
         // CSSCAL
         // Add a test for alpha equal to zero.
         SA = 0.0E0
         for (I = 1; I <= 5; I++) { // 100
            MWPCT(I) = (0.0E0,0.0E0)
            MWPCS(I) = (1.0E0,1.0E0)
         } // 100
         csscal(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to one.
         SA = 1.0E0
         for (I = 1; I <= 5; I++) { // 120
            MWPCT(I) = CX(I)
            MWPCS(I) = CX(I)
         } // 120
         csscal(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
         // Add a test for alpha equal to minus one.
         SA = -1.0E0
         for (I = 1; I <= 5; I++) { // 140
            MWPCT(I) = -CX(I)
            MWPCS(I) = -CX(I)
         } // 140
         csscal(5,SA,CX,INCX);
         ctest(5,CX,MWPCT,MWPCS,SFAC);
      }
      RETURN

      // End of CHECK1

      }
      SUBROUTINE CHECK2(SFAC)
      // .. Parameters ..
      int               NOUT;
      const             NOUT=6;
      // .. Scalar Arguments ..
      REAL              SFAC
      // .. Scalars in Common ..
      int               ICASE, INCX, INCY, MODE, N;
      bool              PASS;
      // .. Local Scalars ..
      COMPLEX           CA
      int               I, J, KI, KN, KSIZE, LENX, LENY, LINCX, LINCY, MX, MY;
      // .. Local Arrays ..
      COMPLEX           CDOT(1), CSIZE1(4), CSIZE2(7,2), CSIZE3(14), CT10X(7,4,4), CT10Y(7,4,4), CT6(4,4), CT7(4,4), CT8(7,4,4), CTY0(1), CX(7), CX0(1), CX1(7), CY(7), CY0(1), CY1(7)
      int               INCXS(4), INCYS(4), LENS(4,2), NS(4);
      // .. External Functions ..
      COMPLEX           CDOTC, CDOTU
      // EXTERNAL CDOTC, CDOTU
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSWAP, CTEST
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Data statements ..
      DATA              CA/(0.4E0,-0.7E0)/
      DATA              INCXS/1, 2, -2, -1/
      DATA              INCYS/1, -2, 1, -2/
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/
      DATA              NS/0, 1, 2, 4/
      DATA              CX1/(0.7E0,-0.8E0), (-0.4E0,-0.7E0), (-0.1E0,-0.9E0), (0.2E0,-0.8E0), (-0.9E0,-0.4E0), (0.1E0,0.4E0), (-0.6E0,0.6E0)/
      DATA              CY1/(0.6E0,-0.6E0), (-0.9E0,0.5E0), (0.7E0,-0.6E0), (0.1E0,-0.5E0), (-0.1E0,-0.2E0), (-0.5E0,-0.3E0), (0.8E0,-0.7E0)/
      DATA              ((CT8(I,J,1),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (-1.55E0,0.5E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (-1.55E0,0.5E0), (0.03E0,-0.89E0), (-0.38E0,-0.96E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT8(I,J,2),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.07E0,-0.89E0), (-0.9E0,0.5E0), (0.42E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.78E0,0.06E0), (-0.9E0,0.5E0), (0.06E0,-0.13E0), (0.1E0,-0.5E0), (-0.77E0,-0.49E0), (-0.5E0,-0.3E0), (0.52E0,-1.51E0)/
      DATA              ((CT8(I,J,3),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.07E0,-0.89E0), (-1.18E0,-0.31E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.78E0,0.06E0), (-1.54E0,0.97E0), (0.03E0,-0.89E0), (-0.18E0,-1.31E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT8(I,J,4),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (-0.9E0,0.5E0), (0.05E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0), (-0.9E0,0.5E0), (0.05E0,-0.6E0), (0.1E0,-0.5E0), (-0.77E0,-0.49E0), (-0.5E0,-0.3E0), (0.32E0,-1.16E0)/
      DATA              CT7/(0.0E0,0.0E0), (-0.06E0,-0.90E0), (0.65E0,-0.47E0), (-0.34E0,-1.22E0), (0.0E0,0.0E0), (-0.06E0,-0.90E0), (-0.59E0,-1.46E0), (-1.04E0,-0.04E0), (0.0E0,0.0E0), (-0.06E0,-0.90E0), (-0.83E0,0.59E0), (0.07E0,-0.37E0), (0.0E0,0.0E0), (-0.06E0,-0.90E0), (-0.76E0,-1.15E0), (-1.33E0,-1.82E0)/
      DATA              CT6/(0.0E0,0.0E0), (0.90E0,0.06E0), (0.91E0,-0.77E0), (1.80E0,-0.10E0), (0.0E0,0.0E0), (0.90E0,0.06E0), (1.45E0,0.74E0), (0.20E0,0.90E0), (0.0E0,0.0E0), (0.90E0,0.06E0), (-0.55E0,0.23E0), (0.83E0,-0.39E0), (0.0E0,0.0E0), (0.90E0,0.06E0), (1.04E0,0.79E0), (1.95E0,1.22E0)/
      DATA              ((CT10X(I,J,1),I=1,7),J=1,4)/(0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (-0.9E0,0.5E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (-0.9E0,0.5E0), (0.7E0,-0.6E0), (0.1E0,-0.5E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10X(I,J,2),I=1,7),J=1,4)/(0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.6E0), (-0.4E0,-0.7E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.8E0,-0.7E0), (-0.4E0,-0.7E0), (-0.1E0,-0.2E0), (0.2E0,-0.8E0), (0.7E0,-0.6E0), (0.1E0,0.4E0), (0.6E0,-0.6E0)/
      DATA              ((CT10X(I,J,3),I=1,7),J=1,4)/(0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.9E0,0.5E0), (-0.4E0,-0.7E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.1E0,-0.5E0), (-0.4E0,-0.7E0), (0.7E0,-0.6E0), (0.2E0,-0.8E0), (-0.9E0,0.5E0), (0.1E0,0.4E0), (0.6E0,-0.6E0)/
      DATA              ((CT10X(I,J,4),I=1,7),J=1,4)/(0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.7E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.7E0,-0.6E0), (-0.1E0,-0.2E0), (0.8E0,-0.7E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,1),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.4E0,-0.7E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.4E0,-0.7E0), (-0.1E0,-0.9E0), (0.2E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,2),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.1E0,-0.9E0), (-0.9E0,0.5E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.6E0,0.6E0), (-0.9E0,0.5E0), (-0.9E0,-0.4E0), (0.1E0,-0.5E0), (-0.1E0,-0.9E0), (-0.5E0,-0.3E0), (0.7E0,-0.8E0)/
      DATA              ((CT10Y(I,J,3),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.1E0,-0.9E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.6E0,0.6E0), (-0.9E0,-0.4E0), (-0.1E0,-0.9E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,4),I=1,7),J=1,4)/(0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.9E0,0.5E0), (-0.4E0,-0.7E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.9E0,0.5E0), (-0.4E0,-0.7E0), (0.1E0,-0.5E0), (-0.1E0,-0.9E0), (-0.5E0,-0.3E0), (0.2E0,-0.8E0)/
      DATA              CSIZE1/(0.0E0,0.0E0), (0.9E0,0.9E0), (1.63E0,1.73E0), (2.90E0,2.78E0)/
      DATA              CSIZE3/(0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (1.17E0,1.17E0), (1.17E0,1.17E0), (1.17E0,1.17E0), (1.17E0,1.17E0), (1.17E0,1.17E0), (1.17E0,1.17E0), (1.17E0,1.17E0)/
      DATA              CSIZE2/(0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0), (1.54E0,1.54E0), (1.54E0,1.54E0), (1.54E0,1.54E0), (1.54E0,1.54E0), (1.54E0,1.54E0), (1.54E0,1.54E0), (1.54E0,1.54E0)/
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
            if (ICASE == 1) {
               // .. CDOTC ..
               CDOT(1) = CDOTC(N,CX,INCX,CY,INCY)
               ctest(1,CDOT,CT6(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE == 2) {
               // .. CDOTU ..
               CDOT(1) = CDOTU(N,CX,INCX,CY,INCY)
               ctest(1,CDOT,CT7(KN,KI),CSIZE1(KN),SFAC);
            } else if (ICASE == 3) {
               // .. CAXPY ..
               caxpy(N,CA,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC);
            } else if (ICASE == 4) {
               // .. CCOPY ..
               ccopy(N,CX,INCX,CY,INCY);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0E0);
               if (KI == 1) {
                  CX0(1) = (42.0E0,43.0E0)
                  CY0(1) = (44.0E0,45.0E0)
                  if (N == 0) {
                     CTY0(1) = CY0(1)
                  } else {
                     CTY0(1) = CX0(1)
                  }
                  LINCX = INCX
                  INCX = 0
                  LINCY = INCY
                  INCY = 0
                  ccopy(N,CX0,INCX,CY0,INCY);
                  ctest(1,CY0,CTY0,CSIZE3,1.0E0);
                  INCX = LINCX
                  INCY = LINCY
               }
            } else if (ICASE == 5) {
               // .. CSWAP ..
               cswap(N,CX,INCX,CY,INCY);
               ctest(LENX,CX,CT10X(1,KN,KI),CSIZE3,1.0E0);
               ctest(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0E0);
            } else {
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            }

         } // 40
      } // 60
      RETURN

      // End of CHECK2

      }
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC)
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
      REAL             SFAC
      int              LEN;
      // .. Array Arguments ..
      REAL             SCOMP(LEN), SSIZE(LEN), STRUE(LEN)
      // .. Scalars in Common ..
      int              ICASE, INCX, INCY, MODE, N;
      bool             PASS;
      // .. Local Scalars ..
      REAL             SD
      int              I;
      // .. External Functions ..
      REAL             SDIFF
      // EXTERNAL SDIFF
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // .. Common blocks ..
      // COMMON /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
      // .. Executable Statements ..

      for (I = 1; I <= LEN; I++) { // 40
         SD = SCOMP(I) - STRUE(I)
         IF (ABS(SFAC*SD) .LE. ABS(SSIZE(I))*EPSILON(ZERO)) GO TO 40

                              // HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).

         if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
         PASS = false;
         WRITE (NOUT,99999)
         WRITE (NOUT,99998)
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, I, SCOMP(I), STRUE(I), SD, SSIZE(I)
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE  I                            ', ' COMP(I)                             TRUE(I)  DIFFERENCE', '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)

      // End of STEST

      }
      SUBROUTINE STEST1(SCOMP1,STRUE1,SSIZE,SFAC)
      // ************************* STEST1 *****************************

      // THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
      // REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
      // ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      REAL              SCOMP1, SFAC, STRUE1
      // .. Array Arguments ..
      REAL              SSIZE(*)
      // .. Local Arrays ..
      REAL              SCOMP(1), STRUE(1)
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Executable Statements ..

      SCOMP(1) = SCOMP1
      STRUE(1) = STRUE1
      stest(1,SCOMP,STRUE,SSIZE,SFAC);

      RETURN

      // End of STEST1

      }
      REAL             FUNCTION SDIFF(SA,SB)
      // ********************************* SDIFF **************************
      // COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15

      // .. Scalar Arguments ..
      REAL                            SA, SB
      // .. Executable Statements ..
      SDIFF = SA - SB
      RETURN

      // End of SDIFF

      }
      SUBROUTINE CTEST(LEN,CCOMP,CTRUE,CSIZE,SFAC)
      // **************************** CTEST *****************************

      // C.L. LAWSON, JPL, 1978 DEC 6

      // .. Scalar Arguments ..
      REAL             SFAC
      int              LEN;
      // .. Array Arguments ..
      COMPLEX          CCOMP(LEN), CSIZE(LEN), CTRUE(LEN)
      // .. Local Scalars ..
      int              I;
      // .. Local Arrays ..
      REAL             SCOMP(20), SSIZE(20), STRUE(20)
      // .. External Subroutines ..
      // EXTERNAL STEST
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
      // .. Executable Statements ..
      for (I = 1; I <= LEN; I++) { // 20
         SCOMP(2*I-1) = REAL(CCOMP(I))
         SCOMP(2*I) = AIMAG(CCOMP(I))
         STRUE(2*I-1) = REAL(CTRUE(I))
         STRUE(2*I) = AIMAG(CTRUE(I))
         SSIZE(2*I-1) = REAL(CSIZE(I))
         SSIZE(2*I) = AIMAG(CSIZE(I))
      } // 20

      stest(2*LEN,SCOMP,STRUE,SSIZE,SFAC);
      RETURN

      // End of CTEST

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
      if (ICOMP == ITRUE) GO TO 40;

                             // HERE ICOMP IS NOT EQUAL TO ITRUE.

      if (.NOT. PASS) GO TO 20;
                              // PRINT FAIL MESSAGE AND HEADER.
      PASS = false;
      WRITE (NOUT,99999)
      WRITE (NOUT,99998)
   20 ID = ICOMP - ITRUE
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, ICOMP, ITRUE, ID
      } // 40
      RETURN

99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE                               ', ' COMP                                TRUE     DIFFERENCE', /1X)
99997 FORMAT (1X,I4,I3,3I5,2I36,I12)

      // End of ITEST1

      }
      SUBROUTINE CB1NRM2(N,INCX,THRESH)
      // Compare NRM2 with a reference computation using combinations
      // of the following values:

      // 0, very small, small, ulp, 1, 1/ulp, big, very big, infinity, NaN

      // one of these values is used to initialize x(1) and x(2:N) is
      // these values.
      // filled with random values from [-1,1] scaled by another of

      // This routine is adapted from the test suite provided by
      // Anderson E. (2017)
      // Algorithm 978: Safe Scaling in the Level 1 BLAS
      // ACM Trans Math Softw 44:1--28
      // https://doi.org/10.1145/3061665

      // .. Scalar Arguments ..
      int               INCX, N;
      REAL              THRESH

*  =====================================================================
      // .. Parameters ..
      int               NMAX, NOUT, NV;
      const             NMAX=20, NOUT=6, NV=10;
      REAL              HALF, ONE, THREE, TWO, ZERO
      const             HALF=0.5E+0, ONE=1.0E+0, TWO= 2.0E+0, THREE=3.0E+0, ZERO=0.0E+0;
      // .. External Functions ..
      REAL              SCNRM2
      // EXTERNAL SCNRM2
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, ABS, CMPLX, MAX, MIN, REAL, SQRT
      // .. Model parameters ..
      REAL              BIGNUM, SAFMAX, SAFMIN, SMLNUM, ULP
      const             BIGNUM=0.1014120480E+32, SAFMAX=0.8507059173E+38, SAFMIN=0.1175494351E-37, SMLNUM=0.9860761315E-31, ULP=0.1192092896E-06;
      // .. Local Scalars ..
      COMPLEX           ROGUE
      REAL              SNRM, TRAT, V0, V1, WORKSSQ, Y1, Y2, YMAX, YMIN, YNRM, ZNRM
      int               I, IV, IW, IX, KS;
      bool              FIRST;
      // .. Local Arrays ..
      COMPLEX           X(NMAX), Z(NMAX)
      REAL              VALUES(NV), WORK(NMAX)
      // .. Executable Statements ..
      VALUES(1) = ZERO
      VALUES(2) = TWO*SAFMIN
      VALUES(3) = SMLNUM
      VALUES(4) = ULP
      VALUES(5) = ONE
      VALUES(6) = ONE / ULP
      VALUES(7) = BIGNUM
      VALUES(8) = SAFMAX
      VALUES(9) = SXVALS(V0,2)
      VALUES(10) = SXVALS(V0,3)
      ROGUE = CMPLX(1234.5678E+0,-1234.5678E+0)
      FIRST = true;

      // Check that the arrays are large enough

      if (N*ABS(INCX).GT.NMAX) {
         WRITE (NOUT,99) "SCNRM2", NMAX, INCX, N, N*ABS(INCX)
         RETURN
      }

      // Zero-sized inputs are tested in STEST1.
      if (N.LE.0) {
         RETURN
      }

      // Generate 2*(N-1) values in (-1,1).

      KS = 2*(N-1)
      for (I = 1; I <= KS; I++) {
         random_number(WORK(I));
         WORK(I) = ONE - TWO*WORK(I)
      }

      // Compute the sum of squares of the random values
      // by an unscaled algorithm.

      WORKSSQ = ZERO
      for (I = 1; I <= KS; I++) {
         WORKSSQ = WORKSSQ + WORK(I)*WORK(I)
      }

      // Construct the test vector with one known value
      // and the rest from the random work array multiplied
      // by a scaling factor.

      for (IV = 1; IV <= NV; IV++) {
         V0 = VALUES(IV)
         if (ABS(V0).GT.ONE) {
            V0 = V0*HALF*HALF
         }
         Z(1) = CMPLX(V0,-THREE*V0)
         for (IW = 1; IW <= NV; IW++) {
            V1 = VALUES(IW)
            if (ABS(V1).GT.ONE) {
               V1 = (V1*HALF) / SQRT(REAL(KS+1))
            }
            for (I = 1; I <= N-1; I++) {
               Z(I+1) = CMPLX(V1*WORK(2*I-1),V1*WORK(2*I))
            }

            // Compute the expected value of the 2-norm

            Y1 = ABS(V0) * SQRT(10.0E0)
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

            if ((Y1 != Y1) || (Y2 != Y2)) {
               // add to propagate NaN
               YNRM = Y1 + Y2
            } else if (YMIN == YMAX) {
               YNRM = SQRT(TWO)*YMAX
            } else if (YMAX == ZERO) {
               YNRM = ZERO
            } else {
               YNRM = YMAX*SQRT(ONE + (YMIN / YMAX)**2)
            }

            // Fill the input array to SCNRM2 with steps of incx

            for (I = 1; I <= N; I++) {
               X(I) = ROGUE
            }
            IX = 1
            if (INCX < 0) IX = 1 - (N-1)*INCX;
            for (I = 1; I <= N; I++) {
               X(IX) = Z(I)
               IX = IX + INCX
            }

            // Call SCNRM2 to compute the 2-norm

            SNRM = SCNRM2(N,X,INCX)

            // Compare SNRM and ZNRM.  Roundoff error grows like O(n)
            // in this implementation so we scale the test ratio accordingly.

            if (INCX == 0) {
               Y1 = ABS(REAL(X(1)))
               Y2 = ABS(AIMAG(X(1)))
               YMIN = MIN(Y1, Y2)
               YMAX = MAX(Y1, Y2)
               if ((Y1 != Y1) || (Y2 != Y2)) {
                  // add to propagate NaN
                  ZNRM = Y1 + Y2
               } else if (YMIN == YMAX) {
                  ZNRM = SQRT(TWO)*YMAX
               } else if (YMAX == ZERO) {
                  ZNRM = ZERO
               } else {
                  ZNRM = YMAX * SQRT(ONE + (YMIN / YMAX)**2)
               }
               ZNRM = SQRT(REAL(n)) * ZNRM
            } else {
               ZNRM = YNRM
            }

            // The tests for NaN rely on the compiler not being overly
            // aggressive and removing the statements altogether.
            if ((SNRM != SNRM) || (ZNRM != ZNRM)) {
               if ((SNRM != SNRM).NEQV.(ZNRM != ZNRM)) {
                  TRAT = ONE / ULP
               } else {
                  TRAT = ZERO
               }
            } else if (ZNRM == ZERO) {
               TRAT = SNRM / ULP
            } else {
               TRAT = (ABS(SNRM-ZNRM) / ZNRM) / (TWO*REAL(N)*ULP)
            }
            if ((TRAT != TRAT) || (TRAT.GE.THRESH)) {
               if (FIRST) {
                  FIRST = false;
                  WRITE(NOUT,99999)
               }
               WRITE (NOUT,98) "SCNRM2", N, INCX, IV, IW, TRAT
            }
         }
      }
99999 FORMAT ('                                       FAIL')
   99 FORMAT ( ' Not enough space to test ', A6, ': NMAX = ',I6, ', INCX = ',I6,/,'   N = ',I6,', must be at least ',I6 )
   98 FORMAT( 1X, A6, ': N=', I6,', INCX=', I4, ', IV=', I2, ', IW=', I2, ', test=', E15.8 )
      RETURN
      CONTAINS
      REAL FUNCTION SXVALS(XX,K)
      // .. Scalar Arguments ..
      REAL              XX
      int               K;
      // .. Local Scalars ..
      REAL              X, Y, YY, Z
      // .. Intrinsic Functions ..
      // INTRINSIC HUGE
      // .. Executable Statements ..
      Y = HUGE(XX)
      Z = YY
      if (K == 1) {
         X = -Z
      } else if (K == 2) {
         X = Z
      } else if (K == 3) {
         X = Z / Z
      }
      SXVALS = X
      RETURN
      }
      }
