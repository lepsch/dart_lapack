      SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA,BETA;
      int     K,LDA,LDC,N;
      String    TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),C(LDC,*)
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
      // INTRINSIC DBLE,DCMPLX,DCONJG,MAX
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
      double           RTEMP;
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      double           ONE,ZERO;
      const     ONE=1.0D+0,ZERO=0.0D+0;
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
      } else if ((.NOT.LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'C'))) {
          INFO = 2
      } else if (N.LT.0) {
          INFO = 3
      } else if (K.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 7
      } else if (LDC.LT.MAX(1,N)) {
          INFO = 10
      }
      if (INFO.NE.0) {
          CALL XERBLA('ZHERK ',INFO)
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
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
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
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              }
          }
          RETURN
      }

      // Start the operations.

      if (LSAME(TRANS,'N')) {

         // Form  C := alpha*A*A**H + beta*C.

          if (UPPER) {
              DO 130 J = 1,N
                  if (BETA.EQ.ZERO) {
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 100 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
                  } else {
                      C(J,J) = DBLE(C(J,J))
                  }
                  DO 120 L = 1,K
                      if (A(J,L).NE.DCMPLX(ZERO)) {
                          TEMP = ALPHA*DCONJG(A(J,L))
                          DO 110 I = 1,J - 1
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(I,L))
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
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  } else {
                      C(J,J) = DBLE(C(J,J))
                  }
                  DO 170 L = 1,K
                      if (A(J,L).NE.DCMPLX(ZERO)) {
                          TEMP = ALPHA*DCONJG(A(J,L))
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(J,L))
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      }
  170             CONTINUE
  180         CONTINUE
          }
      } else {

         // Form  C := alpha*A**H*A + beta*C.

          if (UPPER) {
              DO 220 J = 1,N
                  DO 200 I = 1,J - 1
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  190                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  200             CONTINUE
                  RTEMP = ZERO
                  DO 210 L = 1,K
                      RTEMP = RTEMP + DBLE(DCONJG(A(L,J))*A(L,J))
  210             CONTINUE
                  if (BETA.EQ.ZERO) {
                      C(J,J) = ALPHA*RTEMP
                  } else {
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  }
  220         CONTINUE
          } else {
              DO 260 J = 1,N
                  RTEMP = ZERO
                  DO 230 L = 1,K
                      RTEMP = RTEMP + DBLE(DCONJG(A(L,J))*A(L,J))
  230             CONTINUE
                  if (BETA.EQ.ZERO) {
                      C(J,J) = ALPHA*RTEMP
                  } else {
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  }
                  DO 250 I = J + 1,N
                      TEMP = ZERO
                      DO 240 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  240                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  250             CONTINUE
  260         CONTINUE
          }
      }

      RETURN

      // End of ZHERK

      }
