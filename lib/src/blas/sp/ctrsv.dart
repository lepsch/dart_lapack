      void ctrsv(UPLO,TRANS,DIAG,N, final Matrix<double> A, final int LDA, X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      Complex A(LDA,*),X(*);
      // ..

      Complex ZERO;
      const     ZERO= (0.0,0.0);
      Complex TEMP;
      int     I,INFO,IX,J,JX,KX;
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
      // INTRINSIC CONJG,MAX
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
      } else if (LDA < max(1,N)) {
          INFO = 6;
      } else if (INCX == 0) {
          INFO = 8;
      }
      if (INFO != 0) {
          xerbla('CTRSV ',INFO);
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

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (lsame(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (lsame(UPLO,'U')) {
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 20
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/A(J,J);
                          TEMP = X(J);
                          for (I = J - 1; I >= 1; I--) { // 10
                              X[I] = X(I) - TEMP*A(I,J);
                          } // 10
                      }
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX;
                  for (J = N; J >= 1; J--) { // 40
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/A(J,J);
                          TEMP = X(JX);
                          IX = JX;
                          for (I = J - 1; I >= 1; I--) { // 30
                              IX = IX - INCX;
                              X[IX] = X(IX) - TEMP*A(I,J);
                          } // 30
                      }
                      JX = JX - INCX;
                  } // 40
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/A(J,J);
                          TEMP = X(J);
                          for (I = J + 1; I <= N; I++) { // 50
                              X[I] = X(I) - TEMP*A(I,J);
                          } // 50
                      }
                  } // 60
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/A(J,J);
                          TEMP = X(JX);
                          IX = JX;
                          for (I = J + 1; I <= N; I++) { // 70
                              IX = IX + INCX;
                              X[IX] = X(IX) - TEMP*A(I,J);
                          } // 70
                      }
                      JX = JX + INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (lsame(UPLO,'U')) {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 110
                      TEMP = X(J);
                      if (NOCONJ) {
                          for (I = 1; I <= J - 1; I++) { // 90
                              TEMP = TEMP - A(I,J)*X(I);
                          } // 90
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          for (I = 1; I <= J - 1; I++) { // 100
                              TEMP = TEMP - CONJG(A(I,J))*X(I);
                          } // 100
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X[J] = TEMP;
                  } // 110
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 140
                      IX = KX;
                      TEMP = X(JX);
                      if (NOCONJ) {
                          for (I = 1; I <= J - 1; I++) { // 120
                              TEMP = TEMP - A(I,J)*X(IX);
                              IX = IX + INCX;
                          } // 120
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          for (I = 1; I <= J - 1; I++) { // 130
                              TEMP = TEMP - CONJG(A(I,J))*X(IX);
                              IX = IX + INCX;
                          } // 130
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X[JX] = TEMP;
                      JX = JX + INCX;
                  } // 140
              }
          } else {
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 170
                      TEMP = X(J);
                      if (NOCONJ) {
                          for (I = N; I >= J + 1; I--) { // 150
                              TEMP = TEMP - A(I,J)*X(I);
                          } // 150
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          for (I = N; I >= J + 1; I--) { // 160
                              TEMP = TEMP - CONJG(A(I,J))*X(I);
                          } // 160
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X[J] = TEMP;
                  } // 170
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 200
                      IX = KX;
                      TEMP = X(JX);
                      if (NOCONJ) {
                          for (I = N; I >= J + 1; I--) { // 180
                              TEMP = TEMP - A(I,J)*X(IX);
                              IX = IX - INCX;
                          } // 180
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          for (I = N; I >= J + 1; I--) { // 190
                              TEMP = TEMP - CONJG(A(I,J))*X(IX);
                              IX = IX - INCX;
                          } // 190
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X[JX] = TEMP;
                      JX = JX - INCX;
                  } // 200
              }
          }
      }

      }
