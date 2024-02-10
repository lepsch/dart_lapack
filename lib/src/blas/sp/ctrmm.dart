      void ctrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA, final Matrix<double> A, final int LDA, B,LDB) {

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex ALPHA;
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      Complex A(LDA,*),B(LDB,*);
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX
      Complex TEMP;
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOCONJ,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      Complex ONE;
      const     ONE= (1.0,0.0);
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      // ..

      // Test the input parameters.

      LSIDE = lsame(SIDE,'L');
      if (LSIDE) {
          NROWA = M;
      } else {
          NROWA = N;
      }
      NOCONJ = lsame(TRANSA,'T');
      NOUNIT = lsame(DIAG,'N');
      UPPER = lsame(UPLO,'U');

      INFO = 0;
      if (( !LSIDE) && ( !lsame(SIDE,'R'))) {
          INFO = 1;
      } else if (( !UPPER) && ( !lsame(UPLO,'L'))) {
          INFO = 2;
      } else if (( !lsame(TRANSA,'N')) && ( !lsame(TRANSA,'T')) && ( !lsame(TRANSA,'C'))) {
          INFO = 3;
      } else if (( !lsame(DIAG,'U')) && ( !lsame(DIAG,'N'))) {
          INFO = 4;
      } else if (M < 0) {
          INFO = 5;
      } else if (N < 0) {
          INFO = 6;
      } else if (LDA < max(1,NROWA)) {
          INFO = 9;
      } else if (LDB < max(1,M)) {
          INFO = 11;
      }
      if (INFO != 0) {
          xerbla('CTRMM ',INFO);
          return;
      }

      // Quick return if possible.

      if (M == 0 || N == 0) return;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          for (J = 1; J <= N; J++) { // 20
              for (I = 1; I <= M; I++) { // 10
                  B[I][J] = ZERO;
              } // 10
          } // 20
          return;
      }

      // Start the operations.

      if (LSIDE) {
          if (lsame(TRANSA,'N')) {

            // Form  B := alpha*A*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 50
                      for (K = 1; K <= M; K++) { // 40
                          if (B(K,J) != ZERO) {
                              TEMP = ALPHA*B(K,J);
                              for (I = 1; I <= K - 1; I++) { // 30
                                  B[I][J] = B(I,J) + TEMP*A(I,K);
                              } // 30
                              if (NOUNIT) TEMP = TEMP*A(K,K);
                              B[K][J] = TEMP;
                          }
                      } // 40
                  } // 50
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      for (K = M; K >= 1; K--) { // 70
                          if (B(K,J) != ZERO) {
                              TEMP = ALPHA*B(K,J);
                              B[K][J] = TEMP;
                              if (NOUNIT) B(K,J) = B(K,J)*A(K,K);
                              for (I = K + 1; I <= M; I++) { // 60
                                  B[I][J] = B(I,J) + TEMP*A(I,K);
                              } // 60
                          }
                      } // 70
                  } // 80
              }
          } else {

            // Form  B := alpha*A**T*B   or   B := alpha*A**H*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 120
                      for (I = M; I >= 1; I--) { // 110
                          TEMP = B(I,J);
                          if (NOCONJ) {
                              if (NOUNIT) TEMP = TEMP*A(I,I);
                              for (K = 1; K <= I - 1; K++) { // 90
                                  TEMP = TEMP + A(K,I)*B(K,J);
                              } // 90
                          } else {
                              if (NOUNIT) TEMP = TEMP*CONJG(A(I,I));
                              for (K = 1; K <= I - 1; K++) { // 100
                                  TEMP = TEMP + CONJG(A(K,I))*B(K,J);
                              } // 100
                          }
                          B[I][J] = ALPHA*TEMP;
                      } // 110
                  } // 120
              } else {
                  for (J = 1; J <= N; J++) { // 160
                      for (I = 1; I <= M; I++) { // 150
                          TEMP = B(I,J);
                          if (NOCONJ) {
                              if (NOUNIT) TEMP = TEMP*A(I,I);
                              for (K = I + 1; K <= M; K++) { // 130
                                  TEMP = TEMP + A(K,I)*B(K,J);
                              } // 130
                          } else {
                              if (NOUNIT) TEMP = TEMP*CONJG(A(I,I));
                              for (K = I + 1; K <= M; K++) { // 140
                                  TEMP = TEMP + CONJG(A(K,I))*B(K,J);
                              } // 140
                          }
                          B[I][J] = ALPHA*TEMP;
                      } // 150
                  } // 160
              }
          }
      } else {
          if (lsame(TRANSA,'N')) {

            // Form  B := alpha*B*A.

              if (UPPER) {
                  for (J = N; J >= 1; J--) { // 200
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = 1; I <= M; I++) { // 170
                          B[I][J] = TEMP*B(I,J);
                      } // 170
                      for (K = 1; K <= J - 1; K++) { // 190
                          if (A(K,J) != ZERO) {
                              TEMP = ALPHA*A(K,J);
                              for (I = 1; I <= M; I++) { // 180
                                  B[I][J] = B(I,J) + TEMP*B(I,K);
                              } // 180
                          }
                      } // 190
                  } // 200
              } else {
                  for (J = 1; J <= N; J++) { // 240
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = 1; I <= M; I++) { // 210
                          B[I][J] = TEMP*B(I,J);
                      } // 210
                      for (K = J + 1; K <= N; K++) { // 230
                          if (A(K,J) != ZERO) {
                              TEMP = ALPHA*A(K,J);
                              for (I = 1; I <= M; I++) { // 220
                                  B[I][J] = B(I,J) + TEMP*B(I,K);
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
                          if (A(J,K) != ZERO) {
                              if (NOCONJ) {
                                  TEMP = ALPHA*A(J,K);
                              } else {
                                  TEMP = ALPHA*CONJG(A(J,K));
                              }
                              for (I = 1; I <= M; I++) { // 250
                                  B[I][J] = B(I,J) + TEMP*B(I,K);
                              } // 250
                          }
                      } // 260
                      TEMP = ALPHA;
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = TEMP*A(K,K);
                          } else {
                              TEMP = TEMP*CONJG(A(K,K));
                          }
                      }
                      if (TEMP != ONE) {
                          for (I = 1; I <= M; I++) { // 270
                              B[I][K] = TEMP*B(I,K);
                          } // 270
                      }
                  } // 280
              } else {
                  for (K = N; K >= 1; K--) { // 320
                      for (J = K + 1; J <= N; J++) { // 300
                          if (A(J,K) != ZERO) {
                              if (NOCONJ) {
                                  TEMP = ALPHA*A(J,K);
                              } else {
                                  TEMP = ALPHA*CONJG(A(J,K));
                              }
                              for (I = 1; I <= M; I++) { // 290
                                  B[I][J] = B(I,J) + TEMP*B(I,K);
                              } // 290
                          }
                      } // 300
                      TEMP = ALPHA;
                      if (NOUNIT) {
                          if (NOCONJ) {
                              TEMP = TEMP*A(K,K);
                          } else {
                              TEMP = TEMP*CONJG(A(K,K));
                          }
                      }
                      if (TEMP != ONE) {
                          for (I = 1; I <= M; I++) { // 310
                              B[I][K] = TEMP*B(I,K);
                          } // 310
                      }
                  } // 320
              }
          }
      }

      }
