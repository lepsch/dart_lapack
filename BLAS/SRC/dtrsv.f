      SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),X(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,J,JX,KX;
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
      // INTRINSIC MAX
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) {
          INFO = 2
      } else if (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,N)) {
          INFO = 6
      } else if (INCX == 0) {
          INFO = 8
      }
      if (INFO.NE.0) {
          xerbla('DTRSV ',INFO);
          RETURN
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      NOUNIT = LSAME(DIAG,'N')

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX.LE.0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX.NE.1) {
          KX = 1
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  DO 20 J = N,1,-1
                      if (X(J).NE.ZERO) {
                          if (NOUNIT) X(J) = X(J)/A(J,J);
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
                          } // 10
                      }
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      if (X(JX).NE.ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/A(J,J);
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
                          } // 30
                      }
                      JX = JX - INCX
                  } // 40
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J).NE.ZERO) {
                          if (NOUNIT) X(J) = X(J)/A(J,J);
                          TEMP = X(J)
                          for (I = J + 1; I <= N; I++) { // 50
                              X(I) = X(I) - TEMP*A(I,J)
                          } // 50
                      }
                  } // 60
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX).NE.ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/A(J,J);
                          TEMP = X(JX)
                          IX = JX
                          for (I = J + 1; I <= N; I++) { // 70
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
                          } // 70
                      }
                      JX = JX + INCX
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 100
                      TEMP = X(J)
                      for (I = 1; I <= J - 1; I++) { // 90
                          TEMP = TEMP - A(I,J)*X(I)
                      } // 90
                      if (NOUNIT) TEMP = TEMP/A(J,J);
                      X(J) = TEMP
                  } // 100
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 120
                      TEMP = X(JX)
                      IX = KX
                      for (I = 1; I <= J - 1; I++) { // 110
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX + INCX
                      } // 110
                      if (NOUNIT) TEMP = TEMP/A(J,J);
                      X(JX) = TEMP
                      JX = JX + INCX
                  } // 120
              }
          } else {
              if (INCX == 1) {
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      DO 130 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(I)
                      } // 130
                      if (NOUNIT) TEMP = TEMP/A(J,J);
                      X(J) = TEMP
                  } // 140
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      DO 150 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX - INCX
                      } // 150
                      if (NOUNIT) TEMP = TEMP/A(J,J);
                      X(JX) = TEMP
                      JX = JX - INCX
                  } // 160
              }
          }
      }

      RETURN

      // End of DTRSV

      }
