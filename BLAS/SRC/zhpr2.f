      SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      int     INCX,INCY,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 AP(*),X(*),Y(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP1,TEMP2
      int     I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE,DCONJG
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      } else if (INCY == 0) {
          INFO = 7
      }
      if (INFO != 0) {
          xerbla('ZHPR2 ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) .OR. (ALPHA == ZERO)) RETURN

      // Set up the start points in X and Y if the increments are not both
      // unity.

      if ((INCX != 1) .OR. (INCY != 1)) {
          if (INCX.GT.0) {
              KX = 1
          } else {
              KX = 1 - (N-1)*INCX
          }
          if (INCY.GT.0) {
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
                  if ((X(J) != ZERO) .OR. (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      K = KK
                      for (I = 1; I <= J - 1; I++) { // 10
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
                      } // 10
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  }
                  KK = KK + J
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  if ((X(JX) != ZERO) .OR. (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      IX = KX
                      IY = KY
                      for (K = KK; K <= KK + J - 2; K++) { // 30
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
                      } // 30
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
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
                  if ((X(J) != ZERO) .OR. (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      AP(KK) = DBLE(AP(KK)) + DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                      K = KK + 1
                      for (I = J + 1; I <= N; I++) { // 50
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
                      } // 50
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  }
                  KK = KK + N - J + 1
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  if ((X(JX) != ZERO) .OR. (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      AP(KK) = DBLE(AP(KK)) + DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                      IX = JX
                      IY = JY
                      for (K = KK + 1; K <= KK + N - J; K++) { // 70
                          IX = IX + INCX
                          IY = IY + INCY
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                      } // 70
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  }
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
              } // 80
          }
      }

      RETURN

      // End of ZHPR2

      }
