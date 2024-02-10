      void stbsv(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,K,LDA,N;
      String    DIAG,TRANS,UPLO;
      double A(LDA,*),X(*);
      // ..

      double ZERO;
      const     ZERO=0.0;
      double TEMP;
      int     I,INFO,IX,J,JX,KPLUS1,KX,L;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX,MIN
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame(UPLO,'U') && !lsame(UPLO,'L')) {
          INFO = 1;
      } else if ( !lsame(TRANS,'N') && !lsame(TRANS,'T') && !lsame(TRANS,'C')) {
          INFO = 2;
      } else if ( !lsame(DIAG,'U') && !lsame(DIAG,'N')) {
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
          xerbla('STBSV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOUNIT = lsame(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of A are
      // accessed by sequentially with one pass through A.

      if (lsame(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (lsame(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 20
                      if (X(J) != ZERO) {
                          L = KPLUS1 - J;
                          if (NOUNIT) X(J) = X(J)/A(KPLUS1,J);
                          TEMP = X(J);
                          for (I = J - 1; I >= max(1,J-K); I--) { // 10
                              X[I] = X(I) - TEMP*A(L+I,J);
                          } // 10
                      }
                  } // 20
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 40
                      KX = KX - INCX;
                      if (X(JX) != ZERO) {
                          IX = KX;
                          L = KPLUS1 - J;
                          if (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J);
                          TEMP = X(JX);
                          for (I = J - 1; I >= max(1,J-K); I--) { // 30
                              X[IX] = X(IX) - TEMP*A(L+I,J);
                              IX = IX - INCX;
                          } // 30
                      }
                      JX = JX - INCX;
                  } // 40
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J) != ZERO) {
                          L = 1 - J;
                          if (NOUNIT) X(J) = X(J)/A(1,J);
                          TEMP = X(J);
                          for (I = J + 1; I <= min(N,J+K); I++) { // 50
                              X[I] = X(I) - TEMP*A(L+I,J);
                          } // 50
                      }
                  } // 60
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 80
                      KX = KX + INCX;
                      if (X(JX) != ZERO) {
                          IX = KX;
                          L = 1 - J;
                          if (NOUNIT) X(JX) = X(JX)/A(1,J);
                          TEMP = X(JX);
                          for (I = J + 1; I <= min(N,J+K); I++) { // 70
                              X[IX] = X(IX) - TEMP*A(L+I,J);
                              IX = IX + INCX;
                          } // 70
                      }
                      JX = JX + INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T)*x.

          if (lsame(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 100
                      TEMP = X(J);
                      L = KPLUS1 - J;
                      for (I = max(1,J-K); I <= J - 1; I++) { // 90
                          TEMP = TEMP - A(L+I,J)*X(I);
                      } // 90
                      if (NOUNIT) TEMP = TEMP/A(KPLUS1,J);
                      X[J] = TEMP;
                  } // 100
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 120
                      TEMP = X(JX);
                      IX = KX;
                      L = KPLUS1 - J;
                      for (I = max(1,J-K); I <= J - 1; I++) { // 110
                          TEMP = TEMP - A(L+I,J)*X(IX);
                          IX = IX + INCX;
                      } // 110
                      if (NOUNIT) TEMP = TEMP/A(KPLUS1,J);
                      X[JX] = TEMP;
                      JX = JX + INCX;
                      if (J > K) KX = KX + INCX;
                  } // 120
              }
          } else {
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 140
                      TEMP = X(J);
                      L = 1 - J;
                      for (I = min(N,J+K); I >= J + 1; I--) { // 130
                          TEMP = TEMP - A(L+I,J)*X(I);
                      } // 130
                      if (NOUNIT) TEMP = TEMP/A(1,J);
                      X[J] = TEMP;
                  } // 140
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 160
                      TEMP = X(JX);
                      IX = KX;
                      L = 1 - J;
                      for (I = min(N,J+K); I >= J + 1; I--) { // 150
                          TEMP = TEMP - A(L+I,J)*X(IX);
                          IX = IX - INCX;
                      } // 150
                      if (NOUNIT) TEMP = TEMP/A(1,J);
                      X[JX] = TEMP;
                      JX = JX - INCX;
                      if ((N-J) >= K) KX = KX - INCX;
                  } // 160
              }
          }
      }

      }
