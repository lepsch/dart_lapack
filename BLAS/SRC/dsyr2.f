      SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,INCY,LDA,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),X(*),Y(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP1,TEMP2;
      int     I,INFO,IX,IY,J,JX,JY,KX,KY;
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
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N < 0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      } else if (INCY == 0) {
          INFO = 7
      } else if (LDA < MAX(1,N)) {
          INFO = 9
      }
      if (INFO != 0) {
          xerbla('DSYR2 ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) || (ALPHA == ZERO)) RETURN

      // Set up the start points in X and Y if the increments are not both
      // unity.

      if ((INCX != 1) || (INCY != 1)) {
          if (INCX.GT.0) {
              KX = 1
          } else {
              KX = 1 - (N-1)*INCX
          }
          if (INCY.GT.0) {
              KY = 1
          } else {
              KY = 1 - (N-1)*INCY
          }
          JX = KX
          JY = KY
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      if (LSAME(UPLO,'U')) {

         // Form  A  when A is stored in the upper triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 20
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      for (I = 1; I <= J; I++) { // 10
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
                      } // 10
                  }
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      for (I = 1; I <= J; I++) { // 30
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
                      } // 30
                  }
                  JX = JX + INCX
                  JY = JY + INCY
              } // 40
          }
      } else {

         // Form  A  when A is stored in the lower triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 60
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      for (I = J; I <= N; I++) { // 50
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
                      } // 50
                  }
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      for (I = J; I <= N; I++) { // 70
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
                      } // 70
                  }
                  JX = JX + INCX
                  JY = JY + INCY
              } // 80
          }
      }

      RETURN

      // End of DSYR2

      }
