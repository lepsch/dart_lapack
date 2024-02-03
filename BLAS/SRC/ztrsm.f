      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
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
      // INTRINSIC DCONJG,MAX
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOCONJ,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      COMPLEX*16 ONE
      const     ONE= (1.0D+0,0.0D+0);
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..

      // Test the input parameters.

      LSIDE = LSAME(SIDE,'L')
      if (LSIDE) {
          NROWA = M
      } else {
          NROWA = N
      }
      NOCONJ = LSAME(TRANSA,'T')
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
          xerbla('ZTRSM ',INFO);
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

            // Form  B := alpha*inv( A )*B.

              if (UPPER) {
                  DO 60 J = 1,N
                      if (ALPHA.NE.ONE) {
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      }
                      DO 50 K = M,1,-1
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          }
   50                 CONTINUE
   60             CONTINUE
              } else {
                  DO 100 J = 1,N
                      if (ALPHA.NE.ONE) {
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      }
                      DO 90 K = 1,M
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          }
   90                 CONTINUE
  100             CONTINUE
              }
          } else {

            // Form  B := alpha*inv( A**T )*B
            // or    B := alpha*inv( A**H )*B.

              if (UPPER) {
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          if (NOCONJ) {
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
  110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          } else {
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  120                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          }
                          B(I,J) = TEMP
  130                 CONTINUE
  140             CONTINUE
              } else {
                  DO 180 J = 1,N
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          if (NOCONJ) {
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
  150                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          } else {
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  160                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          }
                          B(I,J) = TEMP
  170                 CONTINUE
  180             CONTINUE
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*inv( A ).

              if (UPPER) {
                  DO 230 J = 1,N
                      if (ALPHA.NE.ONE) {
                          DO 190 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  190                     CONTINUE
                      }
                      DO 210 K = 1,J - 1
                          if (A(K,J).NE.ZERO) {
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200                         CONTINUE
                          }
  210                 CONTINUE
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          DO 220 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  220                     CONTINUE
                      }
  230             CONTINUE
              } else {
                  DO 280 J = N,1,-1
                      if (ALPHA.NE.ONE) {
                          DO 240 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  240                     CONTINUE
                      }
                      DO 260 K = J + 1,N
                          if (A(K,J).NE.ZERO) {
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250                         CONTINUE
                          }
  260                 CONTINUE
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          DO 270 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  270                     CONTINUE
                      }
  280             CONTINUE
              }
          } else {

            // Form  B := alpha*B*inv( A**T )
            // or    B := alpha*B*inv( A**H ).

              if (UPPER) {
                  DO 330 K = N,1,-1
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = ONE/A(K,K)
                          } else {
                              TEMP = ONE/DCONJG(A(K,K))
                          }
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      }
                      DO 310 J = 1,K - 1
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = A(J,K)
                              } else {
                                  TEMP = DCONJG(A(J,K))
                              }
                              DO 300 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  300                         CONTINUE
                          }
  310                 CONTINUE
                      if (ALPHA.NE.ONE) {
                          DO 320 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  320                     CONTINUE
                      }
  330             CONTINUE
              } else {
                  DO 380 K = 1,N
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = ONE/A(K,K)
                          } else {
                              TEMP = ONE/DCONJG(A(K,K))
                          }
                          DO 340 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  340                     CONTINUE
                      }
                      DO 360 J = K + 1,N
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = A(J,K)
                              } else {
                                  TEMP = DCONJG(A(J,K))
                              }
                              DO 350 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  350                         CONTINUE
                          }
  360                 CONTINUE
                      if (ALPHA.NE.ONE) {
                          DO 370 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  370                     CONTINUE
                      }
  380             CONTINUE
              }
          }
      }

      RETURN

      // End of ZTRSM

      }
