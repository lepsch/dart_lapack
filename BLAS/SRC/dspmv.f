      SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA,BETA;
      int     INCX,INCY,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      double           AP(*),X(*),Y(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ONE,ZERO;
      const     ONE=1.0D+0,ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP1,TEMP2;
      int     I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY;
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
      if (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX.EQ.0) {
          INFO = 6
      } else if (INCY.EQ.0) {
          INFO = 9
      }
      if (INFO.NE.0) {
          xerbla('DSPMV ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN

      // Set up the start points in  X  and  Y.

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

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      // First form  y := beta*y.

      if (BETA.NE.ONE) {
          if (INCY.EQ.1) {
              if (BETA.EQ.ZERO) {
                  for (I = 1; I <= N; I++) { // 10
                      Y(I) = ZERO
   10             CONTINUE
              } else {
                  for (I = 1; I <= N; I++) { // 20
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              }
          } else {
              IY = KY
              if (BETA.EQ.ZERO) {
                  for (I = 1; I <= N; I++) { // 30
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              } else {
                  for (I = 1; I <= N; I++) { // 40
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              }
          }
      }
      IF (ALPHA.EQ.ZERO) RETURN
      KK = 1
      if (LSAME(UPLO,'U')) {

         // Form  y  when AP contains the upper triangle.

          if ((INCX.EQ.1) .AND. (INCY.EQ.1)) {
              for (J = 1; J <= N; J++) { // 60
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  K = KK
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  KK = KK + J
   60         CONTINUE
          } else {
              JX = KX
              JY = KY
              for (J = 1; J <= N; J++) { // 80
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 K = KK,KK + J - 2
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   80         CONTINUE
          }
      } else {

         // Form  y  when AP contains the lower triangle.

          if ((INCX.EQ.1) .AND. (INCY.EQ.1)) {
              for (J = 1; J <= N; J++) { // 100
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*AP(KK)
                  K = KK + 1
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
                  KK = KK + (N-J+1)
  100         CONTINUE
          } else {
              JX = KX
              JY = KY
              for (J = 1; J <= N; J++) { // 120
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*AP(KK)
                  IX = JX
                  IY = JY
                  DO 110 K = KK + 1,KK + N - J
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + (N-J+1)
  120         CONTINUE
          }
      }

      RETURN

      // End of DSPMV

      }
