      SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      int     K,LDA,LDC,N
      String    TRANS,UPLO;
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      bool    LSAME;
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DBLE,DCMPLX,DCONJG,MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      DOUBLE PRECISION RTEMP
      int     I,INFO,J,L,NROWA
      bool    UPPER;
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHERK ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*A**H + beta*C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.DCMPLX(ZERO)) THEN
                          TEMP = ALPHA*DCONJG(A(J,L))
                          DO 110 I = 1,J - 1
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(I,L))
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.DCMPLX(ZERO)) THEN
                          TEMP = ALPHA*DCONJG(A(J,L))
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(J,L))
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**H*A + beta*C.
*
          IF (UPPER) THEN
              DO 220 J = 1,N
                  DO 200 I = 1,J - 1
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
                  RTEMP = ZERO
                  DO 210 L = 1,K
                      RTEMP = RTEMP + DBLE(DCONJG(A(L,J))*A(L,J))
  210             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                      C(J,J) = ALPHA*RTEMP
                  ELSE
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  END IF
  220         CONTINUE
          ELSE
              DO 260 J = 1,N
                  RTEMP = ZERO
                  DO 230 L = 1,K
                      RTEMP = RTEMP + DBLE(DCONJG(A(L,J))*A(L,J))
  230             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                      C(J,J) = ALPHA*RTEMP
                  ELSE
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  END IF
                  DO 250 I = J + 1,N
                      TEMP = ZERO
                      DO 240 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  240                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  250             CONTINUE
  260         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZHERK
*
      END
