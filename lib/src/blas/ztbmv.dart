      void ztbmv(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,K,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      Complex A(LDA,*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      Complex TEMP;
      int     I,INFO,IX,J,JX,KPLUS1,KX,L;
      bool    NOCONJ,NOUNIT;
      // ..
      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG,MAX,MIN
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
          xerbla('ZTBMV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOCONJ = lsame(TRANS,'T');
      NOUNIT = lsame(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX   too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (lsame(TRANS,'N')) {

          // Form  x := A*x.

          if (lsame(UPLO,'U')) {
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

         // Form  x := A**T*x  or  x := A**H*x.

          if (lsame(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 110
                      TEMP = X(J);
                      L = KPLUS1 - J;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(KPLUS1,J);
                          for (I = J - 1; I >= max(1,J-K); I--) { // 90
                              TEMP = TEMP + A(L+I,J)*X(I);
                          } // 90
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(KPLUS1,J));
                          for (I = J - 1; I >= max(1,J-K); I--) { // 100
                              TEMP = TEMP + DCONJG(A(L+I,J))*X(I);
                          } // 100
                      }
                      X[J] = TEMP;
                  } // 110
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 140
                      TEMP = X(JX);
                      KX = KX - INCX;
                      IX = KX;
                      L = KPLUS1 - J;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(KPLUS1,J);
                          for (I = J - 1; I >= max(1,J-K); I--) { // 120
                              TEMP = TEMP + A(L+I,J)*X(IX);
                              IX = IX - INCX;
                          } // 120
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(KPLUS1,J));
                          for (I = J - 1; I >= max(1,J-K); I--) { // 130
                              TEMP = TEMP + DCONJG(A(L+I,J))*X(IX);
                              IX = IX - INCX;
                          } // 130
                      }
                      X[JX] = TEMP;
                      JX = JX - INCX;
                  } // 140
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 170
                      TEMP = X(J);
                      L = 1 - J;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(1,J);
                          for (I = J + 1; I <= min(N,J+K); I++) { // 150
                              TEMP = TEMP + A(L+I,J)*X(I);
                          } // 150
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(1,J));
                          for (I = J + 1; I <= min(N,J+K); I++) { // 160
                              TEMP = TEMP + DCONJG(A(L+I,J))*X(I);
                          } // 160
                      }
                      X[J] = TEMP;
                  } // 170
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 200
                      TEMP = X(JX);
                      KX = KX + INCX;
                      IX = KX;
                      L = 1 - J;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(1,J);
                          for (I = J + 1; I <= min(N,J+K); I++) { // 180
                              TEMP = TEMP + A(L+I,J)*X(IX);
                              IX = IX + INCX;
                          } // 180
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(1,J));
                          for (I = J + 1; I <= min(N,J+K); I++) { // 190
                              TEMP = TEMP + DCONJG(A(L+I,J))*X(IX);
                              IX = IX + INCX;
                          } // 190
                      }
                      X[JX] = TEMP;
                      JX = JX + INCX;
                  } // 200
              }
          }
      }

      return;
      }
