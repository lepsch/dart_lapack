      SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX AP(*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,J,JX,K,KK,KX;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,REAL
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
          CALL XERBLA('CHPR  ',INFO)
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (ALPHA.EQ.REAL(ZERO))) RETURN

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
                      TEMP = ALPHA*CONJG(X(J))
                      K = KK
                      DO 10 I = 1,J - 1
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                      AP(KK+J-1) = REAL(AP(KK+J-1)) + REAL(X(J)*TEMP)
                  } else {
                      AP(KK+J-1) = REAL(AP(KK+J-1))
                  }
                  KK = KK + J
   20         CONTINUE
          } else {
              JX = KX
              DO 40 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      IX = KX
                      DO 30 K = KK,KK + J - 2
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      AP(KK+J-1) = REAL(AP(KK+J-1)) + REAL(X(JX)*TEMP)
                  } else {
                      AP(KK+J-1) = REAL(AP(KK+J-1))
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
                      TEMP = ALPHA*CONJG(X(J))
                      AP(KK) = REAL(AP(KK)) + REAL(TEMP*X(J))
                      K = KK + 1
                      DO 50 I = J + 1,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  } else {
                      AP(KK) = REAL(AP(KK))
                  }
                  KK = KK + N - J + 1
   60         CONTINUE
          } else {
              JX = KX
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      AP(KK) = REAL(AP(KK)) + REAL(TEMP*X(JX))
                      IX = JX
                      DO 70 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          AP(K) = AP(K) + X(IX)*TEMP
   70                 CONTINUE
                  } else {
                      AP(KK) = REAL(AP(KK))
                  }
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          }
      }

      RETURN

      // End of CHPR

      }
