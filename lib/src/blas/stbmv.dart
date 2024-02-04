      void stbmv(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,K,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      REAL ZERO;
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      REAL TEMP;
      int     I,INFO,IX,J,JX,KPLUS1,KX,L;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      //- bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX,MIN
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME(UPLO,'U') && !LSAME(UPLO,'L')) {
          INFO = 1;
      } else if ( !LSAME(TRANS,'N') && !LSAME(TRANS,'T') && !LSAME(TRANS,'C')) {
          INFO = 2;
      } else if ( !LSAME(DIAG,'U') && !LSAME(DIAG,'N')) {
          INFO = 3;
      } else if (N < 0) {
          INFO = 4;
      } else if (K < 0) {
          INFO = 5;
      } else if (LDA < (K+1)) {
          INFO = 7;
      } else if (INCX == 0) {
          INFO = 9;
      }
      if (INFO != 0) {
          xerbla('STBMV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOUNIT = LSAME(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX   too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

          // Form  x := A*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 20
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          L = KPLUS1 - J;
                          for (I = max(1,J-K); I <= J - 1; I++) { // 10
                              X[I] = X(I) + TEMP*A(L+I,J);
                          } // 10
                          if (NOUNIT) X(J) = X(J)*A(KPLUS1,J);
                      }
                  } // 20
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 40
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          L = KPLUS1 - J;
                          for (I = max(1,J-K); I <= J - 1; I++) { // 30
                              X[IX] = X(IX) + TEMP*A(L+I,J);
                              IX = IX + INCX;
                          } // 30
                          if (NOUNIT) X(JX) = X(JX)*A(KPLUS1,J);
                      }
                      JX = JX + INCX;
                      if (J > K) KX = KX + INCX;
                  } // 40
              }
          } else {
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 60
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          L = 1 - J;
                          for (I = min(N,J+K); I >= J + 1; I--) { // 50
                              X[I] = X(I) + TEMP*A(L+I,J);
                          } // 50
                          if (NOUNIT) X(J) = X(J)*A(1,J);
                      }
                  } // 60
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 80
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          L = 1 - J;
                          for (I = min(N,J+K); I >= J + 1; I--) { // 70
                              X[IX] = X(IX) + TEMP*A(L+I,J);
                              IX = IX - INCX;
                          } // 70
                          if (NOUNIT) X(JX) = X(JX)*A(1,J);
                      }
                      JX = JX - INCX;
                      if ((N-J) >= K) KX = KX - INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := A**T*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 100
                      TEMP = X(J);
                      L = KPLUS1 - J;
                      if (NOUNIT) TEMP = TEMP*A(KPLUS1,J);
                      for (I = J - 1; I >= max(1,J-K); I--) { // 90
                          TEMP = TEMP + A(L+I,J)*X(I);
                      } // 90
                      X[J] = TEMP;
                  } // 100
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 120
                      TEMP = X(JX);
                      KX = KX - INCX;
                      IX = KX;
                      L = KPLUS1 - J;
                      if (NOUNIT) TEMP = TEMP*A(KPLUS1,J);
                      for (I = J - 1; I >= max(1,J-K); I--) { // 110
                          TEMP = TEMP + A(L+I,J)*X(IX);
                          IX = IX - INCX;
                      } // 110
                      X[JX] = TEMP;
                      JX = JX - INCX;
                  } // 120
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 140
                      TEMP = X(J);
                      L = 1 - J;
                      if (NOUNIT) TEMP = TEMP*A(1,J);
                      for (I = J + 1; I <= min(N,J+K); I++) { // 130
                          TEMP = TEMP + A(L+I,J)*X(I);
                      } // 130
                      X[J] = TEMP;
                  } // 140
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 160
                      TEMP = X(JX);
                      KX = KX + INCX;
                      IX = KX;
                      L = 1 - J;
                      if (NOUNIT) TEMP = TEMP*A(1,J);
                      for (I = J + 1; I <= min(N,J+K); I++) { // 150
                          TEMP = TEMP + A(L+I,J)*X(IX);
                          IX = IX + INCX;
                      } // 150
                      X[JX] = TEMP;
                      JX = JX + INCX;
                  } // 160
              }
          }
      }

      return;
      }
