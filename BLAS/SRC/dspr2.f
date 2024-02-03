      SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,INCY,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      double           AP(*),X(*),Y(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP1,TEMP2;
      int     I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N < 0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      } else if (INCY == 0) {
          INFO = 7
      }
      if (INFO != 0) {
          xerbla('DSPR2 ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) || (ALPHA == ZERO)) RETURN

      // Set up the start points in X and Y if the increments are not both
      // unity.

      if ((INCX != 1) || (INCY != 1)) {
          if (INCX > 0) {
              KX = 1
          } else {
              KX = 1 - (N-1)*INCX
          }
          if (INCY > 0) {
              KY = 1
          } else {
              KY = 1 - (N-1)*INCY
          }
          JX = KX
          JY = KY
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      if (LSAME(UPLO,'U')) {

         // Form  A  when upper triangle is stored in AP.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 20
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      for (I = 1; I <= J; I++) { // 10
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
                      } // 10
                  }
                  KK = KK + J
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      for (K = KK; K <= KK + J - 1; K++) { // 30
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
                      } // 30
                  }
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
              } // 40
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 60
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      for (I = J; I <= N; I++) { // 50
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
                      } // 50
                  }
                  KK = KK + N - J + 1
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      for (K = KK; K <= KK + N - J; K++) { // 70
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
                      } // 70
                  }
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
              } // 80
          }
      }

      RETURN

      // End of DSPR2

      }
