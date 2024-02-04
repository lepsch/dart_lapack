      void sgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC) {

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double ALPHA,BETA;
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
      // ..
      // .. Array Arguments ..
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
      // ..
      // .. Local Scalars ..
      double TEMP;
      int     I,INFO,J,L,NROWA,NROWB;
      bool    NOTA,NOTB;
      // ..
      // .. Parameters ..
      double ONE,ZERO;
      const     ONE=1.0,ZERO=0.0;
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // transposed and set  NROWA and NROWB  as the number of rows of  A
      // and  B  respectively.

      NOTA = lsame(TRANSA,'N');
      NOTB = lsame(TRANSB,'N');
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
      if (( !NOTA) && ( !lsame(TRANSA,'C')) && ( !lsame(TRANSA,'T'))) {
          INFO = 1;
      } else if (( !NOTB) && ( !lsame(TRANSB,'C')) && ( !lsame(TRANSB,'T'))) {
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
          xerbla('SGEMM ',INFO);
          return;
      }

      // Quick return if possible.

      if ((M == 0) || (N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) return;

      // And if  alpha == zero.

      if (ALPHA == ZERO) {
          if (BETA == ZERO) {
              for (J = 1; J <= N; J++) { // 20
                  for (I = 1; I <= M; I++) { // 10
                      C[I,J] = ZERO;
                  } // 10
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  for (I = 1; I <= M; I++) { // 30
                      C[I,J] = BETA*C(I,J);
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
                          C[I,J] = ZERO;
                      } // 50
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 60
                          C[I,J] = BETA*C(I,J);
                      } // 60
                  }
                  for (L = 1; L <= K; L++) { // 80
                      TEMP = ALPHA*B(L,J);
                      for (I = 1; I <= M; I++) { // 70
                          C[I,J] = C(I,J) + TEMP*A(I,L);
                      } // 70
                  } // 80
              } // 90
          } else {

            // Form  C := alpha*A**T*B + beta*C

              for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 100
                          TEMP = TEMP + A(L,I)*B(L,J);
                      } // 100
                      if (BETA == ZERO) {
                          C[I,J] = ALPHA*TEMP;
                      } else {
                          C[I,J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 110
              } // 120
          }
      } else {
          if (NOTA) {

            // Form  C := alpha*A*B**T + beta*C

              for (J = 1; J <= N; J++) { // 170
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 130
                          C[I,J] = ZERO;
                      } // 130
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 140
                          C[I,J] = BETA*C(I,J);
                      } // 140
                  }
                  for (L = 1; L <= K; L++) { // 160
                      TEMP = ALPHA*B(J,L);
                      for (I = 1; I <= M; I++) { // 150
                          C[I,J] = C(I,J) + TEMP*A(I,L);
                      } // 150
                  } // 160
              } // 170
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              for (J = 1; J <= N; J++) { // 200
                  for (I = 1; I <= M; I++) { // 190
                      TEMP = ZERO;
                      for (L = 1; L <= K; L++) { // 180
                          TEMP = TEMP + A(L,I)*B(J,L);
                      } // 180
                      if (BETA == ZERO) {
                          C[I,J] = ALPHA*TEMP;
                      } else {
                          C[I,J] = ALPHA*TEMP + BETA*C(I,J);
                      }
                  } // 190
              } // 200
          }
      }

      return;
      }