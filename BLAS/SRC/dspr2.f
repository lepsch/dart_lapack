      SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,INCY,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      double           AP(*),X(*),Y(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0D+0;
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
          INFO = 5
      } else if (INCY.EQ.0) {
          INFO = 7
      }
      if (INFO.NE.0) {
          xerbla('DSPR2 ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN

      // Set up the start points in X and Y if the increments are not both
      // unity.

      if ((INCX.NE.1) .OR. (INCY.NE.1)) {
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

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      if (LSAME(UPLO,'U')) {

         // Form  A  when upper triangle is stored in AP.

          if ((INCX.EQ.1) .AND. (INCY.EQ.1)) {
              DO 20 J = 1,N
                  if ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   10                 CONTINUE
                  }
                  KK = KK + J
   20         CONTINUE
          } else {
              DO 40 J = 1,N
                  if ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  }
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   40         CONTINUE
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if ((INCX.EQ.1) .AND. (INCY.EQ.1)) {
              DO 60 J = 1,N
                  if ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) {
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   50                 CONTINUE
                  }
                  KK = KK + N - J + 1
   60         CONTINUE
          } else {
              DO 80 J = 1,N
                  if ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) {
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  }
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
   80         CONTINUE
          }
      }

      RETURN

      // End of DSPR2

      }
