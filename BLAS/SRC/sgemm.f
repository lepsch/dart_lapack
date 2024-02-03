      SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA,BETA
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),B(LDB,*),C(LDC,*)
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
      // INTRINSIC MAX
      // ..
      // .. Local Scalars ..
      REAL TEMP
      int     I,INFO,J,L,NROWA,NROWB;
      bool    NOTA,NOTB;
      // ..
      // .. Parameters ..
      REAL ONE,ZERO
      const     ONE=1.0,ZERO=0.0;
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // transposed and set  NROWA and NROWB  as the number of rows of  A
      // and  B  respectively.

      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      if (NOTA) {
          NROWA = M
      } else {
          NROWA = K
      }
      if (NOTB) {
          NROWB = K
      } else {
          NROWB = N
      }

      // Test the input parameters.

      INFO = 0
      if ((.NOT.NOTA) && (.NOT.LSAME(TRANSA,'C')) && (.NOT.LSAME(TRANSA,'T'))) {
          INFO = 1
      } else if ((.NOT.NOTB) && (.NOT.LSAME(TRANSB,'C')) && (.NOT.LSAME(TRANSB,'T'))) {
          INFO = 2
      } else if (M < 0) {
          INFO = 3
      } else if (N < 0) {
          INFO = 4
      } else if (K < 0) {
          INFO = 5
      } else if (LDA < MAX(1,NROWA)) {
          INFO = 8
      } else if (LDB < MAX(1,NROWB)) {
          INFO = 10
      } else if (LDC < MAX(1,M)) {
          INFO = 13
      }
      if (INFO != 0) {
          xerbla('SGEMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M == 0) || (N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) RETURN

      // And if  alpha == zero.

      if (ALPHA == ZERO) {
          if (BETA == ZERO) {
              for (J = 1; J <= N; J++) { // 20
                  for (I = 1; I <= M; I++) { // 10
                      C(I,J) = ZERO
                  } // 10
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  for (I = 1; I <= M; I++) { // 30
                      C(I,J) = BETA*C(I,J)
                  } // 30
              } // 40
          }
          RETURN
      }

      // Start the operations.

      if (NOTB) {
          if (NOTA) {

            // Form  C := alpha*A*B + beta*C.

              for (J = 1; J <= N; J++) { // 90
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 50
                          C(I,J) = ZERO
                      } // 50
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 60
                          C(I,J) = BETA*C(I,J)
                      } // 60
                  }
                  for (L = 1; L <= K; L++) { // 80
                      TEMP = ALPHA*B(L,J)
                      for (I = 1; I <= M; I++) { // 70
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      } // 70
                  } // 80
              } // 90
          } else {

            // Form  C := alpha*A**T*B + beta*C

              for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 100
                          TEMP = TEMP + A(L,I)*B(L,J)
                      } // 100
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
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
                          C(I,J) = ZERO
                      } // 130
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 140
                          C(I,J) = BETA*C(I,J)
                      } // 140
                  }
                  for (L = 1; L <= K; L++) { // 160
                      TEMP = ALPHA*B(J,L)
                      for (I = 1; I <= M; I++) { // 150
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      } // 150
                  } // 160
              } // 170
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              for (J = 1; J <= N; J++) { // 200
                  for (I = 1; I <= M; I++) { // 190
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 180
                          TEMP = TEMP + A(L,I)*B(J,L)
                      } // 180
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 190
              } // 200
          }
      }

      RETURN

      // End of SGEMM

      }
