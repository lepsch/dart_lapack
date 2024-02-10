      void cherk(UPLO,TRANS,N,K,ALPHA, final Matrix<double> A, final int LDA, BETA,C, final int LDC) {

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double ALPHA,BETA;
      int     K,LDA,LDC,N;
      String    TRANS,UPLO;
      Complex A(LDA,*),C(LDC,*);
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
      // INTRINSIC CMPLX,CONJG,MAX,REAL
      Complex TEMP;
      double RTEMP;
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      double ONE,ZERO;
      const     ONE=1.0,ZERO=0.0;
      // ..

      // Test the input parameters.

      if (lsame(TRANS,'N')) {
          NROWA = N;
      } else {
          NROWA = K;
      }
      UPPER = lsame(UPLO,'U');

      INFO = 0;
      if (( !UPPER) && ( !lsame(UPLO,'L'))) {
          INFO = 1;
      } else if (( !lsame(TRANS,'N')) && ( !lsame(TRANS,'C'))) {
          INFO = 2;
      } else if (N < 0) {
          INFO = 3;
      } else if (K < 0) {
          INFO = 4;
      } else if (LDA < max(1,NROWA)) {
          INFO = 7;
      } else if (LDC < max(1,N)) {
          INFO = 10;
      }
      if (INFO != 0) {
          xerbla('CHERK ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) return;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          if (UPPER) {
              if (BETA == ZERO) {
                  for (J = 1; J <= N; J++) { // 20
                      for (I = 1; I <= J; I++) { // 10
                          C[I][J] = ZERO;
                      } // 10
                  } // 20
              } else {
                  for (J = 1; J <= N; J++) { // 40
                      for (I = 1; I <= J - 1; I++) { // 30
                          C[I][J] = BETA*C(I,J);
                      } // 30
                      C[J][J] = BETA*double(C(J,J));
                  } // 40
              }
          } else {
              if (BETA == ZERO) {
                  for (J = 1; J <= N; J++) { // 60
                      for (I = J; I <= N; I++) { // 50
                          C[I][J] = ZERO;
                      } // 50
                  } // 60
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      C[J][J] = BETA*double(C(J,J));
                      for (I = J + 1; I <= N; I++) { // 70
                          C[I][J] = BETA*C(I,J);
                      } // 70
                  } // 80
              }
          }
          return;
      }

      // Start the operations.

      if (lsame(TRANS,'N')) {

         // Form  C := alpha*A*A**H + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 130
                  if (BETA == ZERO) {
                      for (I = 1; I <= J; I++) { // 90
                          C[I][J] = ZERO;
                      } // 90
                  } else if (BETA != ONE) {
                      for (I = 1; I <= J - 1; I++) { // 100
                          C[I][J] = BETA*C(I,J);
                      } // 100
                      C[J][J] = BETA*double(C(J,J));
                  } else {
                      C[J][J] = double(C(J,J));
                  }
                  for (L = 1; L <= K; L++) { // 120
                      if (A(J,L) != CMPLX(ZERO)) {
                          TEMP = ALPHA*CONJG(A(J,L));
                          for (I = 1; I <= J - 1; I++) { // 110
                              C[I][J] = C(I,J) + TEMP*A(I,L);
                          } // 110
                          C[J][J] = REAL(C(J,J)) + double(TEMP*A(I,L));
                      }
                  } // 120
              } // 130
          } else {
              for (J = 1; J <= N; J++) { // 180
                  if (BETA == ZERO) {
                      for (I = J; I <= N; I++) { // 140
                          C[I][J] = ZERO;
                      } // 140
                  } else if (BETA != ONE) {
                      C[J][J] = BETA*double(C(J,J));
                      for (I = J + 1; I <= N; I++) { // 150
                          C[I][J] = BETA*C(I,J);
                      } // 150
                  } else {
                      C[J][J] = double(C(J,J));
                  }
                  for (L = 1; L <= K; L++) { // 170
                      if (A(J,L) != CMPLX(ZERO)) {
                          TEMP = ALPHA*CONJG(A(J,L));
                          C[J][J] = REAL(C(J,J)) + double(TEMP*A(J,L));
                          for (I = J + 1; I <= N; I++) { // 160
                              C[I][J] = C(I,J) + TEMP*A(I,L);
                          } // 160
                      }
                  } // 170
              } // 180
          }
      } else {

         // Form  C := alpha*A**H*A + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 220
                  for (I = 1; I <= J - 1; I++) { // 200
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 190
                          TEMP = TEMP + CONJG(A(L,I))*A(L,J);
                      } // 190
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 200
                  RTEMP = ZERO;
                  for (L = 1; L <= K; L++) { // 210
                      RTEMP = RTEMP + double(CONJG(A(L,J))*A(L,J));
                  } // 210
                  if (BETA == ZERO) {
                      C[J][J] = ALPHA*RTEMP;
                  } else {
                      C[J][J] = ALPHA*RTEMP + BETA*double(C(J,J));
                  }
              } // 220
          } else {
              for (J = 1; J <= N; J++) { // 260
                  RTEMP = ZERO;
                  for (L = 1; L <= K; L++) { // 230
                      RTEMP = RTEMP + double(CONJG(A(L,J))*A(L,J));
                  } // 230
                  if (BETA == ZERO) {
                      C[J][J] = ALPHA*RTEMP;
                  } else {
                      C[J][J] = ALPHA*RTEMP + BETA*double(C(J,J));
                  }
                  for (I = J + 1; I <= N; I++) { // 250
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 240
                          TEMP = TEMP + CONJG(A(L,I))*A(L,J);
                      } // 240
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 250
              } // 260
          }
      }

      }
