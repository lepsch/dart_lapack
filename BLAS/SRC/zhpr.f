      void zhpr(UPLO,N,ALPHA,X,INCX,AP) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      Complex AP(*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      Complex TEMP;
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
      // INTRINSIC DBLE,DCONJG
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
          xerbla('ZHPR  ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (ALPHA == DBLE(ZERO))) return;

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
                      TEMP = ALPHA*DCONJG(X(J));
                      K = KK;
                      for (I = 1; I <= J - 1; I++) { // 10
                          AP(K) = AP(K) + X(I)*TEMP;
                          K = K + 1;
                      } // 10
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP);
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1));
                  }
                  KK = KK + J;
              } // 20
          } else {
              JX = KX;
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX));
                      IX = KX;
                      for (K = KK; K <= KK + J - 2; K++) { // 30
                          AP(K) = AP(K) + X(IX)*TEMP;
                          IX = IX + INCX;
                      } // 30
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP);
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1));
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
                      TEMP = ALPHA*DCONJG(X(J));
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(J));
                      K = KK + 1;
                      for (I = J + 1; I <= N; I++) { // 50
                          AP(K) = AP(K) + X(I)*TEMP;
                          K = K + 1;
                      } // 50
                  } else {
                      AP(KK) = DBLE(AP(KK));
                  }
                  KK = KK + N - J + 1;
              } // 60
          } else {
              JX = KX;
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX));
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(JX));
                      IX = JX;
                      for (K = KK + 1; K <= KK + N - J; K++) { // 70
                          IX = IX + INCX;
                          AP(K) = AP(K) + X(IX)*TEMP;
                      } // 70
                  } else {
                      AP(KK) = DBLE(AP(KK));
                  }
                  JX = JX + INCX;
                  KK = KK + N - J + 1;
              } // 80
          }
      }

      return;
      }
