      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA,BETA;
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
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
      double           TEMP;
      int     I,INFO,J,L,NROWA,NROWB;
      bool    NOTA,NOTB;
      // ..
      // .. Parameters ..
      double           ONE,ZERO;
      const     ONE=1.0D+0,ZERO=0.0D+0;
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // transposed and set  NROWA and NROWB  as the number of rows of  A
      // and  B  respectively.

      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      if (NOTA) {
          NROWA = M
      } else {
          NROWA = K
      }
      if (NOTB) {
          NROWB = K
      } else {
          NROWB = N
      }

      // Test the input parameters.

      INFO = 0
      if ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. (.NOT.LSAME(TRANSA,'T'))) {
          INFO = 1
      } else if ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. (.NOT.LSAME(TRANSB,'T'))) {
          INFO = 2
      } else if (M.LT.0) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (K.LT.0) {
          INFO = 5
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 8
      } else if (LDB.LT.MAX(1,NROWB)) {
          INFO = 10
      } else if (LDC.LT.MAX(1,M)) {
          INFO = 13
      }
      if (INFO.NE.0) {
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN

      // And if  alpha.eq.zero.

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

      if (NOTB) {
          if (NOTA) {

            // Form  C := alpha*A*B + beta*C.

              DO 90 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  }
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          } else {

            // Form  C := alpha*A**T*B + beta*C

              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  110             CONTINUE
  120         CONTINUE
          }
      } else {
          if (NOTA) {

            // Form  C := alpha*A*B**T + beta*C

              DO 170 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  }
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  190             CONTINUE
  200         CONTINUE
          }
      }

      RETURN

      // End of DGEMM

      }
