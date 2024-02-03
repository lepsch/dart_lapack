      SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      REAL AP(*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL ZERO
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      REAL TEMP
      int     I,INFO,IX,J,JX,K,KK,KX;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..

      // Test the input parameters.

      INFO = 0
      if ( !LSAME(UPLO,'U') && !LSAME(UPLO,'L')) {
          INFO = 1
      } else if ( !LSAME(TRANS,'N') && !LSAME(TRANS,'T') && !LSAME(TRANS,'C')) {
          INFO = 2
      } else if ( !LSAME(DIAG,'U') && !LSAME(DIAG,'N')) {
          INFO = 3
      } else if (N < 0) {
          INFO = 4
      } else if (INCX == 0) {
          INFO = 7
      }
      if (INFO != 0) {
          xerbla('STPSV ',INFO);
          RETURN
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      NOUNIT = LSAME(DIAG,'N')

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX != 1) {
          KX = 1
      }

      // Start the operations. In this version the elements of AP are
      // accessed sequentially with one pass through AP.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              KK = (N* (N+1))/2
              if (INCX == 1) {
                  DO 20 J = N,1,-1
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J)
                          K = KK - 1
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*AP(K)
                              K = K - 1
                          } // 10
                      }
                      KK = KK - J
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX)
                          IX = JX
                          DO 30 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*AP(K)
                          } // 30
                      }
                      JX = JX - INCX
                      KK = KK - J
                  } // 40
              }
          } else {
              KK = 1
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J)
                          K = KK + 1
                          for (I = J + 1; I <= N; I++) { // 50
                              X(I) = X(I) - TEMP*AP(K)
                              K = K + 1
                          } // 50
                      }
                      KK = KK + (N-J+1)
                  } // 60
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX)
                          IX = JX
                          for (K = KK + 1; K <= KK + N - J; K++) { // 70
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*AP(K)
                          } // 70
                      }
                      JX = JX + INCX
                      KK = KK + (N-J+1)
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x.

          if (LSAME(UPLO,'U')) {
              KK = 1
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 100
                      TEMP = X(J)
                      K = KK
                      for (I = 1; I <= J - 1; I++) { // 90
                          TEMP = TEMP - AP(K)*X(I)
                          K = K + 1
                      } // 90
                      if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      X(J) = TEMP
                      KK = KK + J
                  } // 100
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 120
                      TEMP = X(JX)
                      IX = KX
                      for (K = KK; K <= KK + J - 2; K++) { // 110
                          TEMP = TEMP - AP(K)*X(IX)
                          IX = IX + INCX
                      } // 110
                      if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + J
                  } // 120
              }
          } else {
              KK = (N* (N+1))/2
              if (INCX == 1) {
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      K = KK
                      DO 130 I = N,J + 1,-1
                          TEMP = TEMP - AP(K)*X(I)
                          K = K - 1
                      } // 130
                      if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      X(J) = TEMP
                      KK = KK - (N-J+1)
                  } // 140
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      DO 150 K = KK,KK - (N- (J+1)),-1
                          TEMP = TEMP - AP(K)*X(IX)
                          IX = IX - INCX
                      } // 150
                      if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - (N-J+1)
                  } // 160
              }
          }
      }

      RETURN

      // End of STPSV

      }
