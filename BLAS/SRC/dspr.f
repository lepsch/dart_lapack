      SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      double           AP(*),X(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO=0.0D+0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,J,JX,K,KK,KX;
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
      }
      if (INFO.NE.0) {
          xerbla('DSPR  ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN

      // Set the start point in X if the increment is not unity.

      if (INCX.LE.0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX.NE.1) {
          KX = 1
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      if (LSAME(UPLO,'U')) {

         // Form  A  when upper triangle is stored in AP.

          if (INCX.EQ.1) {
              DO 20 J = 1,N
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*X(J)
                      K = KK
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                  }
                  KK = KK + J
   20         CONTINUE
          } else {
              JX = KX
              DO 40 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  }
                  JX = JX + INCX
                  KK = KK + J
   40         CONTINUE
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if (INCX.EQ.1) {
              DO 60 J = 1,N
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*X(J)
                      K = KK
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  }
                  KK = KK + N - J + 1
   60         CONTINUE
          } else {
              JX = KX
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  }
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          }
      }

      RETURN

      // End of DSPR

      }
