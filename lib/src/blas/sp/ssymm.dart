      void ssymm(final int SIDE, final int UPLO, final int M, final int N, final int ALPHA, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int BETA, final int C, final int LDC,) {
  final A = A_.dim();
  final B = B_.dim();

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double ALPHA,BETA;
      int     LDA,LDB,LDC,M,N;
      String    SIDE,UPLO;
      double A(LDA,*),B(LDB,*),C(LDC,*);
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
      // INTRINSIC MAX
      double TEMP1,TEMP2;
      int     I,INFO,J,K,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      double ONE,ZERO;
      const     ONE=1.0,ZERO=0.0;
      // ..

      // Set NROWA as the number of rows of A.

      if (lsame(SIDE,'L')) {
          NROWA = M;
      } else {
          NROWA = N;
      }
      UPPER = lsame(UPLO,'U');

      // Test the input parameters.

      INFO = 0;
      if (( !lsame(SIDE,'L')) && ( !lsame(SIDE,'R'))) {
          INFO = 1;
      } else if (( !UPPER) && ( !lsame(UPLO,'L'))) {
          INFO = 2;
      } else if (M < 0) {
          INFO = 3;
      } else if (N < 0) {
          INFO = 4;
      } else if (LDA < max(1,NROWA)) {
          INFO = 7;
      } else if (LDB < max(1,M)) {
          INFO = 9;
      } else if (LDC < max(1,M)) {
          INFO = 12;
      }
      if (INFO != 0) {
          xerbla('SSYMM ',INFO);
          return;
      }

      // Quick return if possible.

      if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          if (BETA == ZERO) {
              for (J = 1; J <= N; J++) { // 20
                  for (I = 1; I <= M; I++) { // 10
                      C[I][J] = ZERO;
                  } // 10
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  for (I = 1; I <= M; I++) { // 30
                      C[I][J] = BETA*C(I,J);
                  } // 30
              } // 40
          }
          return;
      }

      // Start the operations.

      if (lsame(SIDE,'L')) {

         // Form  C := alpha*A*B + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 70
                  for (I = 1; I <= M; I++) { // 60
                      TEMP1 = ALPHA*B(I,J);
                      TEMP2 = ZERO;
                      for (K = 1; K <= I - 1; K++) { // 50
                          C[K][J] = C(K,J) + TEMP1*A(K,I);
                          TEMP2 = TEMP2 + B(K,J)*A(K,I);
                      } // 50
                      if (BETA == ZERO) {
                          C[I][J] = TEMP1*A(I,I) + ALPHA*TEMP2;
                      } else {
                          C[I][J] = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2;
                      }
                  } // 60
              } // 70
          } else {
              for (J = 1; J <= N; J++) { // 100
                  for (I = M; I >= 1; I--) { // 90
                      TEMP1 = ALPHA*B(I,J);
                      TEMP2 = ZERO;
                      for (K = I + 1; K <= M; K++) { // 80
                          C[K][J] = C(K,J) + TEMP1*A(K,I);
                          TEMP2 = TEMP2 + B(K,J)*A(K,I);
                      } // 80
                      if (BETA == ZERO) {
                          C[I][J] = TEMP1*A(I,I) + ALPHA*TEMP2;
                      } else {
                          C[I][J] = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2;
                      }
                  } // 90
              } // 100
          }
      } else {

         // Form  C := alpha*B*A + beta*C.

          for (J = 1; J <= N; J++) { // 170
              TEMP1 = ALPHA*A(J,J);
              if (BETA == ZERO) {
                  for (I = 1; I <= M; I++) { // 110
                      C[I][J] = TEMP1*B(I,J);
                  } // 110
              } else {
                  for (I = 1; I <= M; I++) { // 120
                      C[I][J] = BETA*C(I,J) + TEMP1*B(I,J);
                  } // 120
              }
              for (K = 1; K <= J - 1; K++) { // 140
                  if (UPPER) {
                      TEMP1 = ALPHA*A(K,J);
                  } else {
                      TEMP1 = ALPHA*A(J,K);
                  }
                  for (I = 1; I <= M; I++) { // 130
                      C[I][J] = C(I,J) + TEMP1*B(I,K);
                  } // 130
              } // 140
              for (K = J + 1; K <= N; K++) { // 160
                  if (UPPER) {
                      TEMP1 = ALPHA*A(J,K);
                  } else {
                      TEMP1 = ALPHA*A(K,J);
                  }
                  for (I = 1; I <= M; I++) { // 150
                      C[I][J] = C(I,J) + TEMP1*B(I,K);
                  } // 150
              } // 160
          } // 170
      }

      }
