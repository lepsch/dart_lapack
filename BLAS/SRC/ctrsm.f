      SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*)
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
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOCONJ,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      COMPLEX ONE
      const     ONE= (1.0E+0,0.0E+0);
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
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
          xerbla('CTRSM ',INFO);
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

            // Form  B := alpha*inv( A )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 60
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 30
                              B(I,J) = ALPHA*B(I,J)
                          } // 30
                      }
                      DO 50 K = M,1,-1
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
                              } // 40
                          }
                      } // 50
                  } // 60
              } else {
                  for (J = 1; J <= N; J++) { // 100
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 70
                              B(I,J) = ALPHA*B(I,J)
                          } // 70
                      }
                      for (K = 1; K <= M; K++) { // 90
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
                              } // 80
                          }
                      } // 90
                  } // 100
              }
          } else {

            // Form  B := alpha*inv( A**T )*B
            // or    B := alpha*inv( A**H )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 140
                      for (I = 1; I <= M; I++) { // 130
                          TEMP = ALPHA*B(I,J)
                          if (NOCONJ) {
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
                              } // 110
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          } else {
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - CONJG(A(K,I))*B(K,J)
                              } // 120
                              IF (NOUNIT) TEMP = TEMP/CONJG(A(I,I))
                          }
                          B(I,J) = TEMP
                      } // 130
                  } // 140
              } else {
                  for (J = 1; J <= N; J++) { // 180
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          if (NOCONJ) {
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
                              } // 150
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          } else {
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - CONJG(A(K,I))*B(K,J)
                              } // 160
                              IF (NOUNIT) TEMP = TEMP/CONJG(A(I,I))
                          }
                          B(I,J) = TEMP
                      } // 170
                  } // 180
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*inv( A ).

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 230
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 190
                              B(I,J) = ALPHA*B(I,J)
                          } // 190
                      }
                      DO 210 K = 1,J - 1
                          if (A(K,J).NE.ZERO) {
                              for (I = 1; I <= M; I++) { // 200
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
                              } // 200
                          }
                      } // 210
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          for (I = 1; I <= M; I++) { // 220
                              B(I,J) = TEMP*B(I,J)
                          } // 220
                      }
                  } // 230
              } else {
                  DO 280 J = N,1,-1
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 240
                              B(I,J) = ALPHA*B(I,J)
                          } // 240
                      }
                      DO 260 K = J + 1,N
                          if (A(K,J).NE.ZERO) {
                              for (I = 1; I <= M; I++) { // 250
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
                              } // 250
                          }
                      } // 260
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          for (I = 1; I <= M; I++) { // 270
                              B(I,J) = TEMP*B(I,J)
                          } // 270
                      }
                  } // 280
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
                              TEMP = ONE/CONJG(A(K,K))
                          }
                          for (I = 1; I <= M; I++) { // 290
                              B(I,K) = TEMP*B(I,K)
                          } // 290
                      }
                      DO 310 J = 1,K - 1
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = A(J,K)
                              } else {
                                  TEMP = CONJG(A(J,K))
                              }
                              for (I = 1; I <= M; I++) { // 300
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
                              } // 300
                          }
                      } // 310
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 320
                              B(I,K) = ALPHA*B(I,K)
                          } // 320
                      }
                  } // 330
              } else {
                  for (K = 1; K <= N; K++) { // 380
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = ONE/A(K,K)
                          } else {
                              TEMP = ONE/CONJG(A(K,K))
                          }
                          for (I = 1; I <= M; I++) { // 340
                              B(I,K) = TEMP*B(I,K)
                          } // 340
                      }
                      DO 360 J = K + 1,N
                          if (A(J,K).NE.ZERO) {
                              if (NOCONJ) {
                                  TEMP = A(J,K)
                              } else {
                                  TEMP = CONJG(A(J,K))
                              }
                              for (I = 1; I <= M; I++) { // 350
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
                              } // 350
                          }
                      } // 360
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 370
                              B(I,K) = ALPHA*B(I,K)
                          } // 370
                      }
                  } // 380
              }
          }
      }

      RETURN

      // End of CTRSM

      }
