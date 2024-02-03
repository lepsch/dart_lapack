      SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      int     K,LDA,LDB,LDC,N;
      String    TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
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
      COMPLEX*16 TEMP1,TEMP2
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      COMPLEX*16 ONE
      const     ONE= (1.0D+0,0.0D+0);
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..

      // Test the input parameters.

      if (LSAME(TRANS,'N')) {
          NROWA = N
      } else {
          NROWA = K
      }
      UPPER = LSAME(UPLO,'U')

      INFO = 0
      if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 1
      } else if ((.NOT.LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'T'))) {
          INFO = 2
      } else if (N.LT.0) {
          INFO = 3
      } else if (K.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 7
      } else if (LDB.LT.MAX(1,NROWA)) {
          INFO = 9
      } else if (LDC.LT.MAX(1,N)) {
          INFO = 12
      }
      if (INFO.NE.0) {
          CALL XERBLA('ZSYR2K',INFO)
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          if (UPPER) {
              if (BETA.EQ.ZERO) {
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              } else {
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              }
          } else {
              if (BETA.EQ.ZERO) {
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              } else {
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              }
          }
          RETURN
      }

      // Start the operations.

      if (LSAME(TRANS,'N')) {

         // Form  C := alpha*A*B**T + alpha*B*A**T + C.

          if (UPPER) {
              DO 130 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  }
                  DO 120 L = 1,K
                      if ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) {
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  110                     CONTINUE
                      }
  120             CONTINUE
  130         CONTINUE
          } else {
              DO 180 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  }
                  DO 170 L = 1,K
                      if ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) {
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  160                     CONTINUE
                      }
  170             CONTINUE
  180         CONTINUE
          }
      } else {

         // Form  C := alpha*A**T*B + alpha*B**T*A + C.

          if (UPPER) {
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + ALPHA*TEMP2
                      }
  200             CONTINUE
  210         CONTINUE
          } else {
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + ALPHA*TEMP2
                      }
  230             CONTINUE
  240         CONTINUE
          }
      }

      RETURN

      // End of ZSYR2K

      }
