      SUBROUTINE CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC);

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA,BETA;
      int     LDA,LDB,LDC,M,N;
      String    SIDE,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*),C(LDC,*);
      // ..

// =====================================================================

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
      COMPLEX TEMP1,TEMP2;
      int     I,INFO,J,K,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      COMPLEX ONE;
      const     ONE= (1.0,0.0);
      COMPLEX ZERO;
      const     ZERO= (0.0,0.0);
      // ..

      // Set NROWA as the number of rows of A.

      if (LSAME(SIDE,'L')) {
          NROWA = M;
      } else {
          NROWA = N;
      }
      UPPER = LSAME(UPLO,'U');

      // Test the input parameters.

      INFO = 0;
      if (( !LSAME(SIDE,'L')) && ( !LSAME(SIDE,'R'))) {
          INFO = 1;
      } else if (( !UPPER) && ( !LSAME(UPLO,'L'))) {
          INFO = 2;
      } else if (M < 0) {
          INFO = 3;
      } else if (N < 0) {
          INFO = 4;
      } else if (LDA < MAX(1,NROWA)) {
          INFO = 7;
      } else if (LDB < MAX(1,M)) {
          INFO = 9;
      } else if (LDC < MAX(1,M)) {
          INFO = 12;
      }
      if (INFO != 0) {
          xerbla('CSYMM ',INFO);
          return;
      }

      // Quick return if possible.

      IF ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) RETURN;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          if (BETA == ZERO) {
              for (J = 1; J <= N; J++) { // 20
                  for (I = 1; I <= M; I++) { // 10
                      C(I,J) = ZERO;
                  } // 10
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  for (I = 1; I <= M; I++) { // 30
                      C(I,J) = BETA*C(I,J);
                  } // 30
              } // 40
          }
          return;
      }

      // Start the operations.

      if (LSAME(SIDE,'L')) {

         // Form  C := alpha*A*B + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 70
                  for (I = 1; I <= M; I++) { // 60
                      TEMP1 = ALPHA*B(I,J);
                      TEMP2 = ZERO;
                      for (K = 1; K <= I - 1; K++) { // 50
                          C(K,J) = C(K,J) + TEMP1*A(K,I);
                          TEMP2 = TEMP2 + B(K,J)*A(K,I);
                      } // 50
                      if (BETA == ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2;
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2;
                      }
                  } // 60
              } // 70
          } else {
              for (J = 1; J <= N; J++) { // 100
                  DO 90 I = M,1,-1;
                      TEMP1 = ALPHA*B(I,J);
                      TEMP2 = ZERO;
                      for (K = I + 1; K <= M; K++) { // 80
                          C(K,J) = C(K,J) + TEMP1*A(K,I);
                          TEMP2 = TEMP2 + B(K,J)*A(K,I);
                      } // 80
                      if (BETA == ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2;
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2;
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
                      C(I,J) = TEMP1*B(I,J);
                  } // 110
              } else {
                  for (I = 1; I <= M; I++) { // 120
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J);
                  } // 120
              }
              for (K = 1; K <= J - 1; K++) { // 140
                  if (UPPER) {
                      TEMP1 = ALPHA*A(K,J);
                  } else {
                      TEMP1 = ALPHA*A(J,K);
                  }
                  for (I = 1; I <= M; I++) { // 130
                      C(I,J) = C(I,J) + TEMP1*B(I,K);
                  } // 130
              } // 140
              for (K = J + 1; K <= N; K++) { // 160
                  if (UPPER) {
                      TEMP1 = ALPHA*A(J,K);
                  } else {
                      TEMP1 = ALPHA*A(K,J);
                  }
                  for (I = 1; I <= M; I++) { // 150
                      C(I,J) = C(I,J) + TEMP1*B(I,K);
                  } // 150
              } // 160
          } // 170
      }

      return;

      // End of CSYMM

      }
