      SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA;
      REAL BETA;
      int     K,LDA,LDB,LDC,N;
      String    TRANS,UPLO;
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
      // INTRINSIC CONJG,MAX,REAL
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP1,TEMP2;
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      REAL ONE;
      const     ONE=1.0;
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
      } else if (( !LSAME(TRANS,'N')) && ( !LSAME(TRANS,'C'))) {
          INFO = 2;
      } else if (N < 0) {
          INFO = 3;
      } else if (K < 0) {
          INFO = 4;
      } else if (LDA < MAX(1,NROWA)) {
          INFO = 7;
      } else if (LDB < MAX(1,NROWA)) {
          INFO = 9;
      } else if (LDC < MAX(1,N)) {
          INFO = 12;
      }
      if (INFO != 0) {
          xerbla('CHER2K',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) RETURN;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          if (UPPER) {
              if (BETA == REAL(ZERO)) {
                  for (J = 1; J <= N; J++) { // 20
                      for (I = 1; I <= J; I++) { // 10
                          C(I,J) = ZERO;
                      } // 10
                  } // 20
              } else {
                  for (J = 1; J <= N; J++) { // 40
                      for (I = 1; I <= J - 1; I++) { // 30
                          C(I,J) = BETA*C(I,J);
                      } // 30
                      C(J,J) = BETA*REAL(C(J,J));
                  } // 40
              }
          } else {
              if (BETA == REAL(ZERO)) {
                  for (J = 1; J <= N; J++) { // 60
                      for (I = J; I <= N; I++) { // 50
                          C(I,J) = ZERO;
                      } // 50
                  } // 60
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      C(J,J) = BETA*REAL(C(J,J));
                      for (I = J + 1; I <= N; I++) { // 70
                          C(I,J) = BETA*C(I,J);
                      } // 70
                  } // 80
              }
          }
          return;
      }

      // Start the operations.

      if (LSAME(TRANS,'N')) {

         // Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
                    // C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 130
                  if (BETA == REAL(ZERO)) {
                      for (I = 1; I <= J; I++) { // 90
                          C(I,J) = ZERO;
                      } // 90
                  } else if (BETA != ONE) {
                      for (I = 1; I <= J - 1; I++) { // 100
                          C(I,J) = BETA*C(I,J);
                      } // 100
                      C(J,J) = BETA*REAL(C(J,J));
                  } else {
                      C(J,J) = REAL(C(J,J));
                  }
                  for (L = 1; L <= K; L++) { // 120
                      if ((A(J,L) != ZERO) || (B(J,L) != ZERO)) {
                          TEMP1 = ALPHA*CONJG(B(J,L));
                          TEMP2 = CONJG(ALPHA*A(J,L));
                          for (I = 1; I <= J - 1; I++) { // 110
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2;
                          } // 110
                          C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2);
                      }
                  } // 120
              } // 130
          } else {
              for (J = 1; J <= N; J++) { // 180
                  if (BETA == REAL(ZERO)) {
                      for (I = J; I <= N; I++) { // 140
                          C(I,J) = ZERO;
                      } // 140
                  } else if (BETA != ONE) {
                      for (I = J + 1; I <= N; I++) { // 150
                          C(I,J) = BETA*C(I,J);
                      } // 150
                      C(J,J) = BETA*REAL(C(J,J));
                  } else {
                      C(J,J) = REAL(C(J,J));
                  }
                  for (L = 1; L <= K; L++) { // 170
                      if ((A(J,L) != ZERO) || (B(J,L) != ZERO)) {
                          TEMP1 = ALPHA*CONJG(B(J,L));
                          TEMP2 = CONJG(ALPHA*A(J,L));
                          for (I = J + 1; I <= N; I++) { // 160
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2;
                          } // 160
                          C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2);
                      }
                  } // 170
              } // 180
          }
      } else {

         // Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
                    // C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 210
                  for (I = 1; I <= J; I++) { // 200
                      TEMP1 = ZERO;
                      TEMP2 = ZERO;
                      for (L = 1; L <= K; L++) { // 190
                          TEMP1 = TEMP1 + CONJG(A(L,I))*B(L,J);
                          TEMP2 = TEMP2 + CONJG(B(L,I))*A(L,J);
                      } // 190
                      if (I == J) {
                          if (BETA == REAL(ZERO)) {
                              C(J,J) = REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2);
                          } else {
                              C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2);
                          }
                      } else {
                          if (BETA == REAL(ZERO)) {
                              C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2;
                          } else {
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2;
                          }
                      }
                  } // 200
              } // 210
          } else {
              for (J = 1; J <= N; J++) { // 240
                  for (I = J; I <= N; I++) { // 230
                      TEMP1 = ZERO;
                      TEMP2 = ZERO;
                      for (L = 1; L <= K; L++) { // 220
                          TEMP1 = TEMP1 + CONJG(A(L,I))*B(L,J);
                          TEMP2 = TEMP2 + CONJG(B(L,I))*A(L,J);
                      } // 220
                      if (I == J) {
                          if (BETA == REAL(ZERO)) {
                              C(J,J) = REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2);
                          } else {
                              C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2);
                          }
                      } else {
                          if (BETA == REAL(ZERO)) {
                              C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2;
                          } else {
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2;
                          }
                      }
                  } // 230
              } // 240
          }
      }

      return;

      // End of CHER2K

      }
