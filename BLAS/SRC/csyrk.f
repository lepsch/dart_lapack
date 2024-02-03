      SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC);

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA,BETA;
      int     K,LDA,LDC,N;
      String    TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),C(LDC,*);
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
      COMPLEX TEMP;
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      COMPLEX ONE;
      const     ONE= (1.0,0.0);
      COMPLEX ZERO;
      const     ZERO= (0.0,0.0);
      // ..

      // Test the input parameters.

      if (LSAME(TRANS,'N')) {
          NROWA = N;
      } else {
          NROWA = K;
      }
      UPPER = LSAME(UPLO,'U');

      INFO = 0;
      if (( !UPPER) && ( !LSAME(UPLO,'L'))) {
          INFO = 1;
      } else if (( !LSAME(TRANS,'N')) && ( !LSAME(TRANS,'T'))) {
          INFO = 2;
      } else if (N < 0) {
          INFO = 3;
      } else if (K < 0) {
          INFO = 4;
      } else if (LDA < MAX(1,NROWA)) {
          INFO = 7;
      } else if (LDC < MAX(1,N)) {
          INFO = 10;
      }
      if (INFO != 0) {
          xerbla('CSYRK ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) RETURN;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          if (UPPER) {
              if (BETA == ZERO) {
                  for (J = 1; J <= N; J++) { // 20
                      for (I = 1; I <= J; I++) { // 10
                          C(I,J) = ZERO;
                      } // 10
                  } // 20
              } else {
                  for (J = 1; J <= N; J++) { // 40
                      for (I = 1; I <= J; I++) { // 30
                          C(I,J) = BETA*C(I,J);
                      } // 30
                  } // 40
              }
          } else {
              if (BETA == ZERO) {
                  for (J = 1; J <= N; J++) { // 60
                      for (I = J; I <= N; I++) { // 50
                          C(I,J) = ZERO;
                      } // 50
                  } // 60
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      for (I = J; I <= N; I++) { // 70
                          C(I,J) = BETA*C(I,J);
                      } // 70
                  } // 80
              }
          }
          return;
      }

      // Start the operations.

      if (LSAME(TRANS,'N')) {

         // Form  C := alpha*A*A**T + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 130
                  if (BETA == ZERO) {
                      for (I = 1; I <= J; I++) { // 90
                          C(I,J) = ZERO;
                      } // 90
                  } else if (BETA != ONE) {
                      for (I = 1; I <= J; I++) { // 100
                          C(I,J) = BETA*C(I,J);
                      } // 100
                  }
                  for (L = 1; L <= K; L++) { // 120
                      if (A(J,L) != ZERO) {
                          TEMP = ALPHA*A(J,L);
                          for (I = 1; I <= J; I++) { // 110
                              C(I,J) = C(I,J) + TEMP*A(I,L);
                          } // 110
                      }
                  } // 120
              } // 130
          } else {
              for (J = 1; J <= N; J++) { // 180
                  if (BETA == ZERO) {
                      for (I = J; I <= N; I++) { // 140
                          C(I,J) = ZERO;
                      } // 140
                  } else if (BETA != ONE) {
                      for (I = J; I <= N; I++) { // 150
                          C(I,J) = BETA*C(I,J);
                      } // 150
                  }
                  for (L = 1; L <= K; L++) { // 170
                      if (A(J,L) != ZERO) {
                          TEMP = ALPHA*A(J,L);
                          for (I = J; I <= N; I++) { // 160
                              C(I,J) = C(I,J) + TEMP*A(I,L);
                          } // 160
                      }
                  } // 170
              } // 180
          }
      } else {

         // Form  C := alpha*A**T*A + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 210
                  for (I = 1; I <= J; I++) { // 200
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 190
                          TEMP = TEMP + A(L,I)*A(L,J);
                      } // 190
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP;
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 200
              } // 210
          } else {
              for (J = 1; J <= N; J++) { // 240
                  for (I = J; I <= N; I++) { // 230
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 220
                          TEMP = TEMP + A(L,I)*A(L,J);
                      } // 220
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP;
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 230
              } // 240
          }
      }

      return;

      // End of CSYRK

      }
