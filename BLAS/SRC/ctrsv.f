      void ctrsv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX ZERO;
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP;
      int     I,INFO,IX,J,JX,KX;
      bool    NOCONJ,NOUNIT;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX
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

      NOCONJ = LSAME(TRANS,'T');
      NOUNIT = LSAME(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  DO 20 J = N,1,-1;
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/A(J,J);
                          TEMP = X(J);
                          DO 10 I = J - 1,1,-1;
                              X(I) = X(I) - TEMP*A(I,J);
                          } // 10
                      }
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX;
                  DO 40 J = N,1,-1;
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/A(J,J);
                          TEMP = X(JX);
                          IX = JX;
                          DO 30 I = J - 1,1,-1;
                              IX = IX - INCX;
                              X(IX) = X(IX) - TEMP*A(I,J);
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
                              X(I) = X(I) - TEMP*A(I,J);
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
                              X(IX) = X(IX) - TEMP*A(I,J);
                          } // 70
                      }
                      JX = JX + INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (LSAME(UPLO,'U')) {
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
                      X(J) = TEMP;
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
                      X(JX) = TEMP;
                      JX = JX + INCX;
                  } // 140
              }
          } else {
              if (INCX == 1) {
                  DO 170 J = N,1,-1;
                      TEMP = X(J);
                      if (NOCONJ) {
                          DO 150 I = N,J + 1,-1;
                              TEMP = TEMP - A(I,J)*X(I);
                          } // 150
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          DO 160 I = N,J + 1,-1;
                              TEMP = TEMP - CONJG(A(I,J))*X(I);
                          } // 160
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X(J) = TEMP;
                  } // 170
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  DO 200 J = N,1,-1;
                      IX = KX;
                      TEMP = X(JX);
                      if (NOCONJ) {
                          DO 180 I = N,J + 1,-1;
                              TEMP = TEMP - A(I,J)*X(IX);
                              IX = IX - INCX;
                          } // 180
                          if (NOUNIT) TEMP = TEMP/A(J,J);
                      } else {
                          DO 190 I = N,J + 1,-1;
                              TEMP = TEMP - CONJG(A(I,J))*X(IX);
                              IX = IX - INCX;
                          } // 190
                          if (NOUNIT) TEMP = TEMP/CONJG(A(J,J));
                      }
                      X(JX) = TEMP;
                      JX = JX - INCX;
                  } // 200
              }
          }
      }

      return;
      }
