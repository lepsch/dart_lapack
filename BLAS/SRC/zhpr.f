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
      if (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX.EQ.0) {
          INFO = 5
      }
      if (INFO.NE.0) {
          xerbla('ZHPR  ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (ALPHA.EQ.DBLE(ZERO))) RETURN

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
              for (J = 1; J <= N; J++) { // 20
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      K = KK
                      DO 10 I = 1,J - 1
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  }
                  KK = KK + J
   20         CONTINUE
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 K = KK,KK + J - 2
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP)
                  } else {
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  }
                  JX = JX + INCX
                  KK = KK + J
   40         CONTINUE
          }
      } else {

         // Form  A  when lower triangle is stored in AP.

          if (INCX.EQ.1) {
              for (J = 1; J <= N; J++) { // 60
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(J))
                      K = KK + 1
                      DO 50 I = J + 1,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  }
                  KK = KK + N - J + 1
   60         CONTINUE
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          AP(K) = AP(K) + X(IX)*TEMP
   70                 CONTINUE
                  } else {
                      AP(KK) = DBLE(AP(KK))
                  }
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          }
      }

      RETURN

      // End of ZHPR

      }
