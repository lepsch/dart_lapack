      SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 AP(*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
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
      // INTRINSIC DBLE,DCONJG
      // ..

      // Test the input parameters.

      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHPR  ',INFO)
          RETURN
      END IF

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (ALPHA.EQ.DBLE(ZERO))) RETURN

      // Set the start point in X if the increment is not unity.

      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      IF (LSAME(UPLO,'U')) THEN

         // Form  A  when upper triangle is stored in AP.

          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      K = KK
                      DO 10 I = 1,J - 1
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  KK = KK + J
   20         CONTINUE
          } else {
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 K = KK,KK + J - 2
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  JX = JX + INCX
                  KK = KK + J
   40         CONTINUE
          END IF
      } else {

         // Form  A  when lower triangle is stored in AP.

          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(J))
                      K = KK + 1
                      DO 50 I = J + 1,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  KK = KK + N - J + 1
   60         CONTINUE
          } else {
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          AP(K) = AP(K) + X(IX)*TEMP
   70                 CONTINUE
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          END IF
      END IF

      RETURN

      // End of ZHPR

      }
