      void ctrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) {

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
          xerbla('CTRMV ',INFO);
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

         // Form  x := A*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 20
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          for (I = 1; I <= J - 1; I++) { // 10
                              X(I) = X(I) + TEMP*A(I,J);
                          } // 10
                          if (NOUNIT) X(J) = X(J)*A(J,J);
                      }
                  } // 20
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 40
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          for (I = 1; I <= J - 1; I++) { // 30
                              X(IX) = X(IX) + TEMP*A(I,J);
                              IX = IX + INCX;
                          } // 30
                          if (NOUNIT) X(JX) = X(JX)*A(J,J);
                      }
                      JX = JX + INCX;
                  } // 40
              }
          } else {
              if (INCX == 1) {
                  DO 60 J = N,1,-1;
                      if (X(J) != ZERO) {
                          TEMP = X(J);
                          DO 50 I = N,J + 1,-1;
                              X(I) = X(I) + TEMP*A(I,J);
                          } // 50
                          if (NOUNIT) X(J) = X(J)*A(J,J);
                      }
                  } // 60
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  DO 80 J = N,1,-1;
                      if (X(JX) != ZERO) {
                          TEMP = X(JX);
                          IX = KX;
                          DO 70 I = N,J + 1,-1;
                              X(IX) = X(IX) + TEMP*A(I,J);
                              IX = IX - INCX;
                          } // 70
                          if (NOUNIT) X(JX) = X(JX)*A(J,J);
                      }
                      JX = JX - INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := A**T*x  or  x := A**H*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  DO 110 J = N,1,-1;
                      TEMP = X(J);
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          DO 90 I = J - 1,1,-1;
                              TEMP = TEMP + A(I,J)*X(I);
                          } // 90
                      } else {
                          if (NOUNIT) TEMP = TEMP*CONJG(A(J,J));
                          DO 100 I = J - 1,1,-1;
                              TEMP = TEMP + CONJG(A(I,J))*X(I);
                          } // 100
                      }
                      X(J) = TEMP;
                  } // 110
              } else {
                  JX = KX + (N-1)*INCX;
                  DO 140 J = N,1,-1;
                      TEMP = X(JX);
                      IX = JX;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          DO 120 I = J - 1,1,-1;
                              IX = IX - INCX;
                              TEMP = TEMP + A(I,J)*X(IX);
                          } // 120
                      } else {
                          if (NOUNIT) TEMP = TEMP*CONJG(A(J,J));
                          DO 130 I = J - 1,1,-1;
                              IX = IX - INCX;
                              TEMP = TEMP + CONJG(A(I,J))*X(IX);
                          } // 130
                      }
                      X(JX) = TEMP;
                      JX = JX - INCX;
                  } // 140
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 170
                      TEMP = X(J);
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          for (I = J + 1; I <= N; I++) { // 150
                              TEMP = TEMP + A(I,J)*X(I);
                          } // 150
                      } else {
                          if (NOUNIT) TEMP = TEMP*CONJG(A(J,J));
                          for (I = J + 1; I <= N; I++) { // 160
                              TEMP = TEMP + CONJG(A(I,J))*X(I);
                          } // 160
                      }
                      X(J) = TEMP;
                  } // 170
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 200
                      TEMP = X(JX);
                      IX = JX;
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          for (I = J + 1; I <= N; I++) { // 180
                              IX = IX + INCX;
                              TEMP = TEMP + A(I,J)*X(IX);
                          } // 180
                      } else {
                          if (NOUNIT) TEMP = TEMP*CONJG(A(J,J));
                          for (I = J + 1; I <= N; I++) { // 190
                              IX = IX + INCX;
                              TEMP = TEMP + CONJG(A(I,J))*X(IX);
                          } // 190
                      }
                      X(JX) = TEMP;
                      JX = JX + INCX;
                  } // 200
              }
          }
      }

      return;
      }
