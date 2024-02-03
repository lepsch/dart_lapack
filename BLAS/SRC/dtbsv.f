      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX);

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,K,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),X(*);
      // ..

// =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,J,JX,KPLUS1,KX,L;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      bool    LSAME;
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
          xerbla('DTBSV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      NOUNIT = LSAME(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of A are
      // accessed by sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  DO 20 J = N,1,-1;
                      if (X(J) != ZERO) {
                          L = KPLUS1 - J;
                          if (NOUNIT) X(J) = X(J)/A(KPLUS1,J);
                          TEMP = X(J);
                          DO 10 I = J - 1,MAX(1,J-K),-1;
                              X(I) = X(I) - TEMP*A(L+I,J);
                          } // 10
                      }
                  } // 20
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  DO 40 J = N,1,-1;
                      KX = KX - INCX;
                      if (X(JX) != ZERO) {
                          IX = KX;
                          L = KPLUS1 - J;
                          if (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J);
                          TEMP = X(JX);
                          DO 30 I = J - 1,MAX(1,J-K),-1;
                              X(IX) = X(IX) - TEMP*A(L+I,J);
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
                          DO 50 I = J + 1,MIN(N,J+K);
                              X(I) = X(I) - TEMP*A(L+I,J);
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
                          DO 70 I = J + 1,MIN(N,J+K);
                              X(IX) = X(IX) - TEMP*A(L+I,J);
                              IX = IX + INCX;
                          } // 70
                      }
                      JX = JX + INCX;
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T)*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 100
                      TEMP = X(J);
                      L = KPLUS1 - J;
                      DO 90 I = MAX(1,J-K),J - 1;
                          TEMP = TEMP - A(L+I,J)*X(I);
                      } // 90
                      if (NOUNIT) TEMP = TEMP/A(KPLUS1,J);
                      X(J) = TEMP;
                  } // 100
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 120
                      TEMP = X(JX);
                      IX = KX;
                      L = KPLUS1 - J;
                      DO 110 I = MAX(1,J-K),J - 1;
                          TEMP = TEMP - A(L+I,J)*X(IX);
                          IX = IX + INCX;
                      } // 110
                      if (NOUNIT) TEMP = TEMP/A(KPLUS1,J);
                      X(JX) = TEMP;
                      JX = JX + INCX;
                      if (J > K) KX = KX + INCX;
                  } // 120
              }
          } else {
              if (INCX == 1) {
                  DO 140 J = N,1,-1;
                      TEMP = X(J);
                      L = 1 - J;
                      DO 130 I = MIN(N,J+K),J + 1,-1;
                          TEMP = TEMP - A(L+I,J)*X(I);
                      } // 130
                      if (NOUNIT) TEMP = TEMP/A(1,J);
                      X(J) = TEMP;
                  } // 140
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  DO 160 J = N,1,-1;
                      TEMP = X(JX);
                      IX = KX;
                      L = 1 - J;
                      DO 150 I = MIN(N,J+K),J + 1,-1;
                          TEMP = TEMP - A(L+I,J)*X(IX);
                          IX = IX - INCX;
                      } // 150
                      if (NOUNIT) TEMP = TEMP/A(1,J);
                      X(JX) = TEMP;
                      JX = JX - INCX;
                      if ((N-J) >= K) KX = KX - INCX;
                  } // 160
              }
          }
      }

      return;

      // End of DTBSV

      }
