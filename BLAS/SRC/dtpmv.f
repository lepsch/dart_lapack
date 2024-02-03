      void dtpmv(UPLO,TRANS,DIAG,N,AP,X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      double           AP(*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,J,JX,K,KK,KX;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      //- bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
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
      } else if (INCX == 0) {
          INFO = 7;
      }
      if (INFO != 0) {
          xerbla('DTPMV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOUNIT = LSAME(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of AP are
      // accessed sequentially with one pass through AP.

      if (LSAME(TRANS,'N')) {

         // Form  x:= A*x.

          if (LSAME(UPLO,'U')) {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 20
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          K = KK;
                          for (I = 1; I <= J - 1; I++) { // 10
                              X(I) = X(I) + TEMP*AP(K);
                              K = K + 1;
                          } // 10
                          if (NOUNIT) X(J) = X(J)*AP(KK+J-1);
                      }
                      KK = KK + J;
                  } // 20
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 40
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          for (K = KK; K <= KK + J - 2; K++) { // 30
                              X(IX) = X(IX) + TEMP*AP(K);
                              IX = IX + INCX;
                          } // 30
                          if (NOUNIT) X(JX) = X(JX)*AP(KK+J-1);
                      }
                      JX = JX + INCX;
                      KK = KK + J;
                  } // 40
              }
          } else {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 60
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          K = KK;
                          for (I = N; I >= J + 1; I--) { // 50
                              X(I) = X(I) + TEMP*AP(K);
                              K = K - 1;
                          } // 50
                          if (NOUNIT) X(J) = X(J)*AP(KK-N+J);
                      }
                      KK = KK - (N-J+1);
                  } // 60
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 80
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          for (K = KK; K >= KK - (N- (J+1)); K--) { // 70
                              X(IX) = X(IX) + TEMP*AP(K);
                              IX = IX - INCX;
                          } // 70
                          if (NOUNIT) X(JX) = X(JX)*AP(KK-N+J);
                      }
                      JX = JX - INCX;
                      KK = KK - (N-J+1);
                  } // 80
              }
          }
      } else {

         // Form  x := A**T*x.

          if (LSAME(UPLO,'U')) {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 100
                      TEMP = X(J);
                      if (NOUNIT) TEMP = TEMP*AP(KK);
                      K = KK - 1;
                      for (I = J - 1; I >= 1; I--) { // 90
                          TEMP = TEMP + AP(K)*X(I);
                          K = K - 1;
                      } // 90
                      X(J) = TEMP;
                      KK = KK - J;
                  } // 100
              } else {
                  JX = KX + (N-1)*INCX;
                  for (J = N; J >= 1; J--) { // 120
                      TEMP = X(JX);
                      IX = JX;
                      if (NOUNIT) TEMP = TEMP*AP(KK);
                      for (K = KK - 1; K >= KK - J + 1; K--) { // 110
                          IX = IX - INCX;
                          TEMP = TEMP + AP(K)*X(IX);
                      } // 110
                      X(JX) = TEMP;
                      JX = JX - INCX;
                      KK = KK - J;
                  } // 120
              }
          } else {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 140
                      TEMP = X(J);
                      if (NOUNIT) TEMP = TEMP*AP(KK);
                      K = KK + 1;
                      for (I = J + 1; I <= N; I++) { // 130
                          TEMP = TEMP + AP(K)*X(I);
                          K = K + 1;
                      } // 130
                      X(J) = TEMP;
                      KK = KK + (N-J+1);
                  } // 140
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 160
                      TEMP = X(JX);
                      IX = JX;
                      if (NOUNIT) TEMP = TEMP*AP(KK);
                      for (K = KK + 1; K <= KK + N - J; K++) { // 150
                          IX = IX + INCX;
                          TEMP = TEMP + AP(K)*X(IX);
                      } // 150
                      X(JX) = TEMP;
                      JX = JX + INCX;
                      KK = KK + (N-J+1);
                  } // 160
              }
          }
      }

      return;
      }
