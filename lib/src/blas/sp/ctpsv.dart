      void ctpsv(final int UPLO, final int TRANS, final int DIAG, final int N, final int AP, final int X, final int INCX,) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,N;
      String    DIAG,TRANS,UPLO;
      Complex AP(*),X(*);
      // ..

      Complex ZERO;
      const     ZERO= (0.0,0.0);
      Complex TEMP;
      int     I,INFO,IX,J,JX,K,KK,KX;
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
      // INTRINSIC CONJG
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
      } else if (INCX == 0) {
          INFO = 7;
      }
      if (INFO != 0) {
          xerbla('CTPSV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOCONJ = lsame(TRANS,'T');
      NOUNIT = lsame(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of AP are
      // accessed sequentially with one pass through AP.

      if (lsame(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (lsame(UPLO,'U')) {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 20
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J);
                          K = KK - 1;
                          for (I = J - 1; I >= 1; I--) { // 10
                              X[I] = X(I) - TEMP*AP(K);
                              K = K - 1;
                          } // 10
                      }
                      KK = KK - J;
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX;
                  for (J = N; J >= 1; J--) { // 40
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX);
                          IX = JX;
                          for (K = KK - 1; K >= KK - J + 1; K--) { // 30
                              IX = IX - INCX;
                              X[IX] = X(IX) - TEMP*AP(K);
                          } // 30
                      }
                      JX = JX - INCX;
                      KK = KK - J;
                  } // 40
              }
          } else {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J);
                          K = KK + 1;
                          for (I = J + 1; I <= N; I++) { // 50
                              X[I] = X(I) - TEMP*AP(K);
                              K = K + 1;
                          } // 50
                      }
                      KK = KK + (N-J+1);
                  } // 60
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX);
                          IX = JX;
                          for (K = KK + 1; K <= KK + N - J; K++) { // 70
                              IX = IX + INCX;
                              X[IX] = X(IX) - TEMP*AP(K);
                          } // 70
                      }
                      JX = JX + INCX;
                      KK = KK + (N-J+1);
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (lsame(UPLO,'U')) {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 110
                      TEMP = X(J);
                      K = KK;
                      if (NOCONJ) {
                          for (I = 1; I <= J - 1; I++) { // 90
                              TEMP = TEMP - AP(K)*X(I);
                              K = K + 1;
                          } // 90
                          if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      } else {
                          for (I = 1; I <= J - 1; I++) { // 100
                              TEMP = TEMP - CONJG(AP(K))*X(I);
                              K = K + 1;
                          } // 100
                          if (NOUNIT) TEMP = TEMP/CONJG(AP(KK+J-1));
                      }
                      X[J] = TEMP;
                      KK = KK + J;
                  } // 110
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 140
                      TEMP = X(JX);
                      IX = KX;
                      if (NOCONJ) {
                          for (K = KK; K <= KK + J - 2; K++) { // 120
                              TEMP = TEMP - AP(K)*X(IX);
                              IX = IX + INCX;
                          } // 120
                          if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      } else {
                          for (K = KK; K <= KK + J - 2; K++) { // 130
                              TEMP = TEMP - CONJG(AP(K))*X(IX);
                              IX = IX + INCX;
                          } // 130
                          if (NOUNIT) TEMP = TEMP/CONJG(AP(KK+J-1));
                      }
                      X[JX] = TEMP;
                      JX = JX + INCX;
                      KK = KK + J;
                  } // 140
              }
          } else {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 170
                      TEMP = X(J);
                      K = KK;
                      if (NOCONJ) {
                          for (I = N; I >= J + 1; I--) { // 150
                              TEMP = TEMP - AP(K)*X(I);
                              K = K - 1;
                          } // 150
                          if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      } else {
                          for (I = N; I >= J + 1; I--) { // 160
                              TEMP = TEMP - CONJG(AP(K))*X(I);
                              K = K - 1;
                          } // 160
                          if (NOUNIT) TEMP = TEMP/CONJG(AP(KK-N+J));
                      }
                      X[J] = TEMP;
                      KK = KK - (N-J+1);
                  } // 170
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 200
                      TEMP = X(JX);
                      IX = KX;
                      if (NOCONJ) {
                          for (K = KK; K >= KK - (N- (J+1)); K--) { // 180
                              TEMP = TEMP - AP(K)*X(IX);
                              IX = IX - INCX;
                          } // 180
                          if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      } else {
                          for (K = KK; K >= KK - (N- (J+1)); K--) { // 190
                              TEMP = TEMP - CONJG(AP(K))*X(IX);
                              IX = IX - INCX;
                          } // 190
                          if (NOUNIT) TEMP = TEMP/CONJG(AP(KK-N+J));
                      }
                      X[JX] = TEMP;
                      JX = JX - INCX;
                      KK = KK - (N-J+1);
                  } // 200
              }
          }
      }

      }
