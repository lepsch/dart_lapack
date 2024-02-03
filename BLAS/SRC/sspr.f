      void sspr(UPLO,N,ALPHA,X,INCX,AP) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA;
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      REAL AP(*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      REAL ZERO;
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      REAL TEMP;
      int     I,INFO,IX,J,JX,K,KK,KX;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME(UPLO,'U') && !LSAME(UPLO,'L')) {
          INFO = 1;
      } else if (N < 0) {
          INFO = 2;
      } else if (INCX == 0) {
          INFO = 5;
      }
      if (INFO != 0) {
          xerbla('SSPR  ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (ALPHA == ZERO)) return;

      // Set the start point in X if the increment is not unity.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1;
      if (LSAME(UPLO,'U')) {

         // Form  A  when upper triangle is stored in AP.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 20
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*X(J);
                      K = KK;
                      for (I = 1; I <= J; I++) { // 10
                          AP(K) = AP(K) + X(I)*TEMP;
                          K = K + 1;
                      } // 10
                  }
                  KK = KK + J;
              } // 20
          } else {
              JX = KX;
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*X(JX);
                      IX = KX;
                      for (K = KK; K <= KK + J - 1; K++) { // 30
                          AP(K) = AP(K) + X(IX)*TEMP;
                          IX = IX + INCX;
                      } // 30
                  }
                  JX = JX + INCX;
                  KK = KK + J;
              } // 40
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 60
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*X(J);
                      K = KK;
                      for (I = J; I <= N; I++) { // 50
                          AP(K) = AP(K) + X(I)*TEMP;
                          K = K + 1;
                      } // 50
                  }
                  KK = KK + N - J + 1;
              } // 60
          } else {
              JX = KX;
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*X(JX);
                      IX = JX;
                      for (K = KK; K <= KK + N - J; K++) { // 70
                          AP(K) = AP(K) + X(IX)*TEMP;
                          IX = IX + INCX;
                      } // 70
                  }
                  JX = JX + INCX;
                  KK = KK + N - J + 1;
              } // 80
          }
      }

      return;

      // End of SSPR

      }
