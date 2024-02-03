      SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX AP(*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,J,JX,K,KK,KX;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,REAL
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      }
      if (INFO != 0) {
          xerbla('CHPR  ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) .OR. (ALPHA == REAL(ZERO))) RETURN

      // Set the start point in X if the increment is not unity.

      if (INCX.LE.0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX != 1) {
          KX = 1
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      if (LSAME(UPLO,'U')) {

         // Form  A  when upper triangle is stored in AP.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 20
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*CONJG(X(J))
                      K = KK
                      for (I = 1; I <= J - 1; I++) { // 10
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
                      } // 10
                      AP(KK+J-1) = REAL(AP(KK+J-1)) + REAL(X(J)*TEMP)
                  } else {
                      AP(KK+J-1) = REAL(AP(KK+J-1))
                  }
                  KK = KK + J
              } // 20
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      IX = KX
                      for (K = KK; K <= KK + J - 2; K++) { // 30
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
                      } // 30
                      AP(KK+J-1) = REAL(AP(KK+J-1)) + REAL(X(JX)*TEMP)
                  } else {
                      AP(KK+J-1) = REAL(AP(KK+J-1))
                  }
                  JX = JX + INCX
                  KK = KK + J
              } // 40
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 60
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*CONJG(X(J))
                      AP(KK) = REAL(AP(KK)) + REAL(TEMP*X(J))
                      K = KK + 1
                      for (I = J + 1; I <= N; I++) { // 50
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
                      } // 50
                  } else {
                      AP(KK) = REAL(AP(KK))
                  }
                  KK = KK + N - J + 1
              } // 60
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      AP(KK) = REAL(AP(KK)) + REAL(TEMP*X(JX))
                      IX = JX
                      for (K = KK + 1; K <= KK + N - J; K++) { // 70
                          IX = IX + INCX
                          AP(K) = AP(K) + X(IX)*TEMP
                      } // 70
                  } else {
                      AP(KK) = REAL(AP(KK))
                  }
                  JX = JX + INCX
                  KK = KK + N - J + 1
              } // 80
          }
      }

      RETURN

      // End of CHPR

      }
