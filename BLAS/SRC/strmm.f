      SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),B(LDB,*)
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
      REAL TEMP
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      REAL ONE,ZERO
      const     ONE=1.0E+0,ZERO=0.0E+0;
      // ..

      // Test the input parameters.

      LSIDE = LSAME(SIDE,'L')
      if (LSIDE) {
          NROWA = M
      } else {
          NROWA = N
      }
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')

      INFO = 0
      if ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) {
          INFO = 1
      } else if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 2
      } else if ((.NOT.LSAME(TRANSA,'N')) .AND. (.NOT.LSAME(TRANSA,'T')) .AND. (.NOT.LSAME(TRANSA,'C'))) {
          INFO = 3
      } else if ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) {
          INFO = 4
      } else if (M.LT.0) {
          INFO = 5
      } else if (N.LT.0) {
          INFO = 6
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 9
      } else if (LDB.LT.MAX(1,M)) {
          INFO = 11
      }
      if (INFO.NE.0) {
          xerbla('STRMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF (M.EQ.0 .OR. N.EQ.0) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      }

      // Start the operations.

      if (LSIDE) {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*A*B.

              if (UPPER) {
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          if (B(K,J).NE.ZERO) {
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          }
   40                 CONTINUE
   50             CONTINUE
              } else {
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          if (B(K,J).NE.ZERO) {
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          }
   70                 CONTINUE
   80             CONTINUE
              }
          } else {

            // Form  B := alpha*A**T*B.

              if (UPPER) {
                  DO 110 J = 1,N
                      DO 100 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 90 K = 1,I - 1
                              TEMP = TEMP + A(K,I)*B(K,J)
   90                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  100                 CONTINUE
  110             CONTINUE
              } else {
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 120 K = I + 1,M
                              TEMP = TEMP + A(K,I)*B(K,J)
  120                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  130                 CONTINUE
  140             CONTINUE
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*A.

              if (UPPER) {
                  DO 180 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  150                 CONTINUE
                      DO 170 K = 1,J - 1
                          if (A(K,J).NE.ZERO) {
                              TEMP = ALPHA*A(K,J)
                              DO 160 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  160                         CONTINUE
                          }
  170                 CONTINUE
  180             CONTINUE
              } else {
                  DO 220 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 190 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  190                 CONTINUE
                      DO 210 K = J + 1,N
                          if (A(K,J).NE.ZERO) {
                              TEMP = ALPHA*A(K,J)
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  200                         CONTINUE
                          }
  210                 CONTINUE
  220             CONTINUE
              }
          } else {

            // Form  B := alpha*B*A**T.

              if (UPPER) {
                  DO 260 K = 1,N
                      DO 240 J = 1,K - 1
                          if (A(J,K).NE.ZERO) {
                              TEMP = ALPHA*A(J,K)
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  230                         CONTINUE
                          }
  240                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      if (TEMP.NE.ONE) {
                          DO 250 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  250                     CONTINUE
                      }
  260             CONTINUE
              } else {
                  DO 300 K = N,1,-1
                      DO 280 J = K + 1,N
                          if (A(J,K).NE.ZERO) {
                              TEMP = ALPHA*A(J,K)
                              DO 270 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  270                         CONTINUE
                          }
  280                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      if (TEMP.NE.ONE) {
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      }
  300             CONTINUE
              }
          }
      }

      RETURN

      // End of STRMM

      }
