      SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,LDA,N;
      String    UPLO;
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
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE,DCONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N.LT.0) {
          INFO = 2
      } else if (INCX.EQ.0) {
          INFO = 5
      } else if (LDA.LT.MAX(1,N)) {
          INFO = 7
      }
      if (INFO.NE.0) {
          CALL XERBLA('ZHER  ',INFO)
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

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      if (LSAME(UPLO,'U')) {

         // Form  A  when A is stored in upper triangle.

          if (INCX.EQ.1) {
              DO 20 J = 1,N
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      DO 10 I = 1,J - 1
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(J)*TEMP)
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
   20         CONTINUE
          } else {
              JX = KX
              DO 40 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 I = 1,J - 1
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(JX)*TEMP)
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
                  JX = JX + INCX
   40         CONTINUE
          }
      } else {

         // Form  A  when A is stored in lower triangle.

          if (INCX.EQ.1) {
              DO 60 J = 1,N
                  if (X(J).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(J))
                      DO 50 I = J + 1,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
   60         CONTINUE
          } else {
              JX = KX
              DO 80 J = 1,N
                  if (X(JX).NE.ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 I = J + 1,N
                          IX = IX + INCX
                          A(I,J) = A(I,J) + X(IX)*TEMP
   70                 CONTINUE
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
                  JX = JX + INCX
   80         CONTINUE
          }
      }

      RETURN

      // End of ZHER

      }
