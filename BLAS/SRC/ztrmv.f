      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
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
      // INTRINSIC DCONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (.NOT.LSAME(TRANS,'N') && .NOT.LSAME(TRANS,'T') && .NOT.LSAME(TRANS,'C')) {
          INFO = 2
      } else if (.NOT.LSAME(DIAG,'U') && .NOT.LSAME(DIAG,'N')) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,N)) {
          INFO = 6
      } else if (INCX == 0) {
          INFO = 8
      }
      if (INFO != 0) {
          xerbla('ZTRMV ',INFO);
          RETURN
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX.LE.0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX != 1) {
          KX = 1
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := A*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 20
                      if (X(J) != ZERO) {
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
                      if (X(JX) != ZERO) {
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
              if (INCX == 1) {
                  DO 60 J = N,1,-1
                      if (X(J) != ZERO) {
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
                      if (X(JX) != ZERO) {
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

         // Form  x := A**T*x  or  x := A**H*x.

          if (LSAME(UPLO,'U')) {
              if (INCX == 1) {
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
                          } // 90
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(J,J));
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
                          } // 100
                      }
                      X(J) = TEMP
                  } // 110
              } else {
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
                          } // 120
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(J,J));
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          } // 130
                      }
                      X(JX) = TEMP
                      JX = JX - INCX
                  } // 140
              }
          } else {
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 170
                      TEMP = X(J)
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          for (I = J + 1; I <= N; I++) { // 150
                              TEMP = TEMP + A(I,J)*X(I)
                          } // 150
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(J,J));
                          for (I = J + 1; I <= N; I++) { // 160
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
                          } // 160
                      }
                      X(J) = TEMP
                  } // 170
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 200
                      TEMP = X(JX)
                      IX = JX
                      if (NOCONJ) {
                          if (NOUNIT) TEMP = TEMP*A(J,J);
                          for (I = J + 1; I <= N; I++) { // 180
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
                          } // 180
                      } else {
                          if (NOUNIT) TEMP = TEMP*DCONJG(A(J,J));
                          for (I = J + 1; I <= N; I++) { // 190
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          } // 190
                      }
                      X(JX) = TEMP
                      JX = JX + INCX
                  } // 200
              }
          }
      }

      RETURN

      // End of ZTRMV

      }
