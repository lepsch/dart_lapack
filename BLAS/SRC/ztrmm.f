      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

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
          xerbla('ZTRMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF (M.EQ.0 .OR. N.EQ.0) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          for (J = 1; J <= N; J++) { // 20
              for (I = 1; I <= M; I++) { // 10
                  B(I,J) = ZERO
              } // 10
          } // 20
          RETURN
      }

      // Start the operations.

      if (LSIDE) {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*A*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 50
                      for (K = 1; K <= M; K++) { // 40
                          if (B(K,J).NE.ZERO) {
                              TEMP = ALPHA*B(K,J)
                              for (I = 1; I <= K - 1; I++) { // 30
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
                              } // 30
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          }
                      } // 40
                  } // 50
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      DO 70 K = M,1,-1
                          if (B(K,J).NE.ZERO) {
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              for (I = K + 1; I <= M; I++) { // 60
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
                              } // 60
                          }
                      } // 70
                  } // 80
              }
          } else {

            // Form  B := alpha*A**T*B   or   B := alpha*A**H*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 120
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          if (NOCONJ) {
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              for (K = 1; K <= I - 1; K++) { // 90
                                  TEMP = TEMP + A(K,I)*B(K,J)
                              } // 90
                          } else {
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              for (K = 1; K <= I - 1; K++) { // 100
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
                              } // 100
                          }
                          B(I,J) = ALPHA*TEMP
                      } // 110
                  } // 120
              } else {
                  for (J = 1; J <= N; J++) { // 160
                      for (I = 1; I <= M; I++) { // 150
                          TEMP = B(I,J)
                          if (NOCONJ) {
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              for (K = I + 1; K <= M; K++) { // 130
                                  TEMP = TEMP + A(K,I)*B(K,J)
                              } // 130
                          } else {
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              for (K = I + 1; K <= M; K++) { // 140
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
                              } // 140
                          }
                          B(I,J) = ALPHA*TEMP
                      } // 150
                  } // 160
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*A.

              if (UPPER) {
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      for (I = 1; I <= M; I++) { // 170
                          B(I,J) = TEMP*B(I,J)
                      } // 170
                      for (K = 1; K <= J - 1; K++) { // 190
                          if (A(K,J).NE.ZERO) {
                              TEMP = ALPHA*A(K,J)
                              for (I = 1; I <= M; I++) { // 180
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
                              } // 180
                          }
                      } // 190
                  } // 200
              } else {
                  for (J = 1; J <= N; J++) { // 240
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      for (I = 1; I <= M; I++) { // 210
                          B(I,J) = TEMP*B(I,J)
                      } // 210
                      for (K = J + 1; K <= N; K++) { // 230
                          if (A(K,J).NE.ZERO) {
                              TEMP = ALPHA*A(K,J)
                              for (I = 1; I <= M; I++) { // 220
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
                              } // 220
                          }
                      } // 230
                  } // 240
              }
          } else {

            // Form  B := alpha*B*A**T   or   B := alpha*B*A**H.

              if (UPPER) {
                  for (K = 1; K <= N; K++) { // 280
                      for (J = 1; J <= K - 1; J++) { // 260
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = ALPHA*A(J,K)
                              } else {
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              }
                              for (I = 1; I <= M; I++) { // 250
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
                              } // 250
                          }
                      } // 260
                      TEMP = ALPHA
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = TEMP*A(K,K)
                          } else {
                              TEMP = TEMP*DCONJG(A(K,K))
                          }
                      }
                      if (TEMP.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 270
                              B(I,K) = TEMP*B(I,K)
                          } // 270
                      }
                  } // 280
              } else {
                  DO 320 K = N,1,-1
                      for (J = K + 1; J <= N; J++) { // 300
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = ALPHA*A(J,K)
                              } else {
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              }
                              for (I = 1; I <= M; I++) { // 290
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
                              } // 290
                          }
                      } // 300
                      TEMP = ALPHA
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = TEMP*A(K,K)
                          } else {
                              TEMP = TEMP*DCONJG(A(K,K))
                          }
                      }
                      if (TEMP.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 310
                              B(I,K) = TEMP*B(I,K)
                          } // 310
                      }
                  } // 320
              }
          }
      }

      RETURN

      // End of ZTRMM

      }
