      void cgemm(TRANSA,TRANSB,M,N,K,ALPHA, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, BETA,C, final int LDC) {

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex ALPHA,BETA;
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
      Complex A(LDA,*),B(LDB,*),C(LDC,*);
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
      int     I,INFO,J,L,NROWA,NROWB;
      bool    CONJA,CONJB,NOTA,NOTB;
      // ..
      // .. Parameters ..
      Complex ONE;
      const     ONE= (1.0,0.0);
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // conjugated or transposed, set  CONJA and CONJB  as true if  A  and
      // B  respectively are to be  transposed but  not conjugated  and set
      // NROWA and  NROWB  as the number of rows of  A  and  B  respectively.

      NOTA = lsame(TRANSA,'N');
      NOTB = lsame(TRANSB,'N');
      CONJA = lsame(TRANSA,'C');
      CONJB = lsame(TRANSB,'C');
      if (NOTA) {
          NROWA = M;
      } else {
          NROWA = K;
      }
      if (NOTB) {
          NROWB = K;
      } else {
          NROWB = N;
      }

      // Test the input parameters.

      INFO = 0;
      if (( !NOTA) && ( !CONJA) && ( !lsame(TRANSA,'T'))) {
          INFO = 1;
      } else if (( !NOTB) && ( !CONJB) && ( !lsame(TRANSB,'T'))) {
          INFO = 2;
      } else if (M < 0) {
          INFO = 3;
      } else if (N < 0) {
          INFO = 4;
      } else if (K < 0) {
          INFO = 5;
      } else if (LDA < max(1,NROWA)) {
          INFO = 8;
      } else if (LDB < max(1,NROWB)) {
          INFO = 10;
      } else if (LDC < max(1,M)) {
          INFO = 13;
      }
      if (INFO != 0) {
          xerbla('CGEMM ',INFO);
          return;
      }

      // Quick return if possible.

      if ((M == 0) || (N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) return;

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

      if (NOTB) {
          if (NOTA) {

            // Form  C := alpha*A*B + beta*C.

              for (J = 1; J <= N; J++) { // 90
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 50
                          C[I][J] = ZERO;
                      } // 50
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 60
                          C[I][J] = BETA*C(I,J);
                      } // 60
                  }
                  for (L = 1; L <= K; L++) { // 80
                      TEMP = ALPHA*B(L,J);
                      for (I = 1; I <= M; I++) { // 70
                          C[I][J] = C(I,J) + TEMP*A(I,L);
                      } // 70
                  } // 80
              } // 90
          } else if (CONJA) {

            // Form  C := alpha*A**H*B + beta*C.

              for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 100
                          TEMP = TEMP + CONJG(A(L,I))*B(L,J);
                      } // 100
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 110
              } // 120
          } else {

            // Form  C := alpha*A**T*B + beta*C

              for (J = 1; J <= N; J++) { // 150
                  for (I = 1; I <= M; I++) { // 140
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 130
                          TEMP = TEMP + A(L,I)*B(L,J);
                      } // 130
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 140
              } // 150
          }
      } else if (NOTA) {
          if (CONJB) {

            // Form  C := alpha*A*B**H + beta*C.

              for (J = 1; J <= N; J++) { // 200
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 160
                          C[I][J] = ZERO;
                      } // 160
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 170
                          C[I][J] = BETA*C(I,J);
                      } // 170
                  }
                  for (L = 1; L <= K; L++) { // 190
                      TEMP = ALPHA*CONJG(B(J,L));
                      for (I = 1; I <= M; I++) { // 180
                          C[I][J] = C(I,J) + TEMP*A(I,L);
                      } // 180
                  } // 190
              } // 200
          } else {

            // Form  C := alpha*A*B**T + beta*C

              for (J = 1; J <= N; J++) { // 250
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 210
                          C[I][J] = ZERO;
                      } // 210
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 220
                          C[I][J] = BETA*C(I,J);
                      } // 220
                  }
                  for (L = 1; L <= K; L++) { // 240
                      TEMP = ALPHA*B(J,L);
                      for (I = 1; I <= M; I++) { // 230
                          C[I][J] = C(I,J) + TEMP*A(I,L);
                      } // 230
                  } // 240
              } // 250
          }
      } else if (CONJA) {
          if (CONJB) {

            // Form  C := alpha*A**H*B**H + beta*C.

              for (J = 1; J <= N; J++) { // 280
                  for (I = 1; I <= M; I++) { // 270
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 260
                          TEMP = TEMP + CONJG(A(L,I))*CONJG(B(J,L));
                      } // 260
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 270
              } // 280
          } else {

            // Form  C := alpha*A**H*B**T + beta*C

              for (J = 1; J <= N; J++) { // 310
                  for (I = 1; I <= M; I++) { // 300
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 290
                          TEMP = TEMP + CONJG(A(L,I))*B(J,L);
                      } // 290
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 300
              } // 310
          }
      } else {
          if (CONJB) {

            // Form  C := alpha*A**T*B**H + beta*C

              for (J = 1; J <= N; J++) { // 340
                  for (I = 1; I <= M; I++) { // 330
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 320
                          TEMP = TEMP + A(L,I)*CONJG(B(J,L));
                      } // 320
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 330
              } // 340
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              for (J = 1; J <= N; J++) { // 370
                  for (I = 1; I <= M; I++) { // 360
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 350
                          TEMP = TEMP + A(L,I)*B(J,L);
                      } // 350
                      if (BETA == ZERO) {
                          C[I][J] = ALPHA*TEMP;
                      } else {
                          C[I][J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 360
              } // 370
          }
      }

      }
