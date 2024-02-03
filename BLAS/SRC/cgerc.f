      SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA
      int     INCX,INCY,LDA,M,N;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),X(*),Y(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,J,JY,KX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0
      if (M.LT.0) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX.EQ.0) {
          INFO = 5
      } else if (INCY.EQ.0) {
          INFO = 7
      } else if (LDA.LT.MAX(1,M)) {
          INFO = 9
      }
      if (INFO.NE.0) {
          xerbla('CGERC ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (INCY.GT.0) {
          JY = 1
      } else {
          JY = 1 - (N-1)*INCY
      }
      if (INCX.EQ.1) {
          for (J = 1; J <= N; J++) { // 20
              if (Y(JY).NE.ZERO) {
                  TEMP = ALPHA*CONJG(Y(JY))
                  for (I = 1; I <= M; I++) { // 10
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              }
              JY = JY + INCY
   20     CONTINUE
      } else {
          if (INCX.GT.0) {
              KX = 1
          } else {
              KX = 1 - (M-1)*INCX
          }
          for (J = 1; J <= N; J++) { // 40
              if (Y(JY).NE.ZERO) {
                  TEMP = ALPHA*CONJG(Y(JY))
                  IX = KX
                  for (I = 1; I <= M; I++) { // 30
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              }
              JY = JY + INCY
   40     CONTINUE
      }

      RETURN

      // End of CGERC

      }
