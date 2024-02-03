      SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA,BETA
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
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
      // INTRINSIC CONJG,MAX
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,J,L,NROWA,NROWB;
      bool    CONJA,CONJB,NOTA,NOTB;
      // ..
      // .. Parameters ..
      COMPLEX ONE
      const     ONE= (1.0E+0,0.0E+0);
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // conjugated or transposed, set  CONJA and CONJB  as true if  A  and
      // B  respectively are to be  transposed but  not conjugated  and set
      // NROWA and  NROWB  as the number of rows of  A  and  B  respectively.

      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
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
      if ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND. (.NOT.LSAME(TRANSA,'T'))) {
          INFO = 1
      } else if ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND. (.NOT.LSAME(TRANSB,'T'))) {
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
          CALL XERBLA('CGEMM ',INFO)
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN

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
          } else if (CONJA) {

            // Form  C := alpha*A**H*B + beta*C.

              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + CONJG(A(L,I))*B(L,J)
  100                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  110             CONTINUE
  120         CONTINUE
          } else {

            // Form  C := alpha*A**T*B + beta*C

              DO 150 J = 1,N
                  DO 140 I = 1,M
                      TEMP = ZERO
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  140             CONTINUE
  150         CONTINUE
          }
      } else if (NOTA) {
          if (CONJB) {

            // Form  C := alpha*A*B**H + beta*C.

              DO 200 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  }
                  DO 190 L = 1,K
                      TEMP = ALPHA*CONJG(B(J,L))
                      DO 180 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          } else {

            // Form  C := alpha*A*B**T + beta*C

              DO 250 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  }
                  DO 240 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 230 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  230                 CONTINUE
  240             CONTINUE
  250         CONTINUE
          }
      } else if (CONJA) {
          if (CONJB) {

            // Form  C := alpha*A**H*B**H + beta*C.

              DO 280 J = 1,N
                  DO 270 I = 1,M
                      TEMP = ZERO
                      DO 260 L = 1,K
                          TEMP = TEMP + CONJG(A(L,I))*CONJG(B(J,L))
  260                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  270             CONTINUE
  280         CONTINUE
          } else {

            // Form  C := alpha*A**H*B**T + beta*C

              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + CONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  300             CONTINUE
  310         CONTINUE
          }
      } else {
          if (CONJB) {

            // Form  C := alpha*A**T*B**H + beta*C

              DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*CONJG(B(J,L))
  320                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  330             CONTINUE
  340         CONTINUE
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              DO 370 J = 1,N
                  DO 360 I = 1,M
                      TEMP = ZERO
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  360             CONTINUE
  370         CONTINUE
          }
      }

      RETURN

      // End of CGEMM

      }
