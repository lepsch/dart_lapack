      SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,LDA,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
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
*
      // Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHER  ',INFO)
          RETURN
      END IF
*
      // Quick return if possible.
*
      IF ((N.EQ.0) .OR. (ALPHA.EQ.DBLE(ZERO))) RETURN
*
      // Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
         // Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      DO 10 I = 1,J - 1
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(J)*TEMP)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 I = 1,J - 1
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(JX)*TEMP)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
*
         // Form  A  when A is stored in lower triangle.
*
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(J))
                      DO 50 I = J + 1,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 I = J + 1,N
                          IX = IX + INCX
                          A(I,J) = A(I,J) + X(IX)*TEMP
   70                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
      // End of ZHER
*
      END
