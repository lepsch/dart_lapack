      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA,BETA;
      int     INCX,INCY,LDA,M,N;
      String    TRANS;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),X(*),Y(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ONE,ZERO;
      const     ONE=1.0D+0,ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY;
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
      if (.NOT.LSAME(TRANS,'N') && .NOT.LSAME(TRANS,'T') && .NOT.LSAME(TRANS,'C')) {
          INFO = 1
      } else if (M.LT.0) {
          INFO = 2
      } else if (N.LT.0) {
          INFO = 3
      } else if (LDA.LT.MAX(1,M)) {
          INFO = 6
      } else if (INCX == 0) {
          INFO = 8
      } else if (INCY == 0) {
          INFO = 11
      }
      if (INFO != 0) {
          xerbla('DGEMV ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M == 0) .OR. (N == 0) .OR. ((ALPHA == ZERO) && (BETA == ONE))) RETURN

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if (LSAME(TRANS,'N')) {
          LENX = N
          LENY = M
      } else {
          LENX = M
          LENY = N
      }
      if (INCX.GT.0) {
          KX = 1
      } else {
          KX = 1 - (LENX-1)*INCX
      }
      if (INCY.GT.0) {
          KY = 1
      } else {
          KY = 1 - (LENY-1)*INCY
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      // First form  y := beta*y.

      if (BETA != ONE) {
          if (INCY == 1) {
              if (BETA == ZERO) {
                  for (I = 1; I <= LENY; I++) { // 10
                      Y(I) = ZERO
                  } // 10
              } else {
                  for (I = 1; I <= LENY; I++) { // 20
                      Y(I) = BETA*Y(I)
                  } // 20
              }
          } else {
              IY = KY
              if (BETA == ZERO) {
                  for (I = 1; I <= LENY; I++) { // 30
                      Y(IY) = ZERO
                      IY = IY + INCY
                  } // 30
              } else {
                  for (I = 1; I <= LENY; I++) { // 40
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
                  } // 40
              }
          }
      }
      if (ALPHA == ZERO) RETURN;
      if (LSAME(TRANS,'N')) {

         // Form  y := alpha*A*x + y.

          JX = KX
          if (INCY == 1) {
              for (J = 1; J <= N; J++) { // 60
                  TEMP = ALPHA*X(JX)
                  for (I = 1; I <= M; I++) { // 50
                      Y(I) = Y(I) + TEMP*A(I,J)
                  } // 50
                  JX = JX + INCX
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  for (I = 1; I <= M; I++) { // 70
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
                  } // 70
                  JX = JX + INCX
              } // 80
          }
      } else {

         // Form  y := alpha*A**T*x + y.

          JY = KY
          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 100
                  TEMP = ZERO
                  for (I = 1; I <= M; I++) { // 90
                      TEMP = TEMP + A(I,J)*X(I)
                  } // 90
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
              } // 100
          } else {
              for (J = 1; J <= N; J++) { // 120
                  TEMP = ZERO
                  IX = KX
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
                  } // 110
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
              } // 120
          }
      }

      RETURN

      // End of DGEMV

      }
