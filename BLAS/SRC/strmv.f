      SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL ZERO
      const     ZERO=0.0E+0;
      // ..
      // .. Local Scalars ..
      REAL TEMP
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
      } else if (INCX.EQ.0) {
          INFO = 8
      }
      if (INFO.NE.0) {
          xerbla('STRMV ',INFO);
          RETURN
      }

      // Quick return if possible.

      if (N.EQ.0) RETURN;

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

         // Form  x := A*x.

          if (LSAME(UPLO,'U')) {
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 20
                      if (X(J).NE.ZERO) {
                          TEMP = X(J)
                          for (I = 1; I <= J - 1; I++) { // 10
                              X(I) = X(I) + TEMP*A(I,J)
                          } // 10
                          if (NOUNIT) X(J) = X(J)*A(J,J);
                      }
                  } // 20
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 40
                      if (X(JX).NE.ZERO) {
                          TEMP = X(JX)
                          IX = KX
                          for (I = 1; I <= J - 1; I++) { // 30
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
                          } // 30
                          if (NOUNIT) X(JX) = X(JX)*A(J,J);
                      }
                      JX = JX + INCX
                  } // 40
              }
          } else {
              if (INCX.EQ.1) {
                  DO 60 J = N,1,-1
                      if (X(J).NE.ZERO) {
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
                          } // 50
                          if (NOUNIT) X(J) = X(J)*A(J,J);
                      }
                  } // 60
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      if (X(JX).NE.ZERO) {
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
                          } // 70
                          if (NOUNIT) X(JX) = X(JX)*A(J,J);
                      }
                      JX = JX - INCX
                  } // 80
              }
          }
      } else {

         // Form  x := A**T*x.

          if (LSAME(UPLO,'U')) {
              if (INCX.EQ.1) {
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + A(I,J)*X(I)
                      } // 90
                      X(J) = TEMP
                  } // 100
              } else {
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      DO 110 I = J - 1,1,-1
                          IX = IX - INCX
                          TEMP = TEMP + A(I,J)*X(IX)
                      } // 110
                      X(JX) = TEMP
                      JX = JX - INCX
                  } // 120
              }
          } else {
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 140
                      TEMP = X(J)
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = J + 1; I <= N; I++) { // 130
                          TEMP = TEMP + A(I,J)*X(I)
                      } // 130
                      X(J) = TEMP
                  } // 140
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 160
                      TEMP = X(JX)
                      IX = JX
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = J + 1; I <= N; I++) { // 150
                          IX = IX + INCX
                          TEMP = TEMP + A(I,J)*X(IX)
                      } // 150
                      X(JX) = TEMP
                      JX = JX + INCX
                  } // 160
              }
          }
      }

      RETURN

      // End of STRMV

      }
