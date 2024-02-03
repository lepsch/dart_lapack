      SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA,BETA;
      int     LDA,LDB,LDC,M,N;
      String    SIDE,UPLO;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),B(LDB,*),C(LDC,*);
      // ..

*  =====================================================================

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
      // .. Local Scalars ..
      double           TEMP1,TEMP2;
      int     I,INFO,J,K,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      double           ONE,ZERO;
      const     ONE=1.0D+0,ZERO=0.0D+0;
      // ..

      // Set NROWA as the number of rows of A.

      if (LSAME(SIDE,'L')) {
          NROWA = M
      } else {
          NROWA = N
      }
      UPPER = LSAME(UPLO,'U')

      // Test the input parameters.

      INFO = 0
      if ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) {
          INFO = 1
      } else if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 2
      } else if (M.LT.0) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 7
      } else if (LDB.LT.MAX(1,M)) {
          INFO = 9
      } else if (LDC.LT.MAX(1,M)) {
          INFO = 12
      }
      if (INFO.NE.0) {
          CALL XERBLA('DSYMM ',INFO)
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          if (BETA.EQ.ZERO) {
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          } else {
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          }
          RETURN
      }

      // Start the operations.

      if (LSAME(SIDE,'L')) {

         // Form  C := alpha*A*B + beta*C.

          if (UPPER) {
              DO 70 J = 1,N
                  DO 60 I = 1,M
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   50                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                      }
   60             CONTINUE
   70         CONTINUE
          } else {
              DO 100 J = 1,N
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   80                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                      }
   90             CONTINUE
  100         CONTINUE
          }
      } else {

         // Form  C := alpha*B*A + beta*C.

          DO 170 J = 1,N
              TEMP1 = ALPHA*A(J,J)
              if (BETA.EQ.ZERO) {
                  DO 110 I = 1,M
                      C(I,J) = TEMP1*B(I,J)
  110             CONTINUE
              } else {
                  DO 120 I = 1,M
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  120             CONTINUE
              }
              DO 140 K = 1,J - 1
                  if (UPPER) {
                      TEMP1 = ALPHA*A(K,J)
                  } else {
                      TEMP1 = ALPHA*A(J,K)
                  }
                  DO 130 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  130             CONTINUE
  140         CONTINUE
              DO 160 K = J + 1,N
                  if (UPPER) {
                      TEMP1 = ALPHA*A(J,K)
                  } else {
                      TEMP1 = ALPHA*A(K,J)
                  }
                  DO 150 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
      }

      RETURN

      // End of DSYMM

      }
