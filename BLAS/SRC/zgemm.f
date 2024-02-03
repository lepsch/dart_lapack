      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      int     K,LDA,LDB,LDC,M,N;
      String    TRANSA,TRANSB;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
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
      // INTRINSIC DCONJG,MAX
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
      int     I,INFO,J,L,NROWA,NROWB;
      bool    CONJA,CONJB,NOTA,NOTB;
      // ..
      // .. Parameters ..
      COMPLEX*16 ONE
      const     ONE= (1.0D+0,0.0D+0);
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..

      // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
      // conjugated or transposed, set  CONJA and CONJB  as true if  A  and
      // B  respectively are to be  transposed but  not conjugated  and set
      // NROWA and NROWB  as the number of rows  of  A  and  B  respectively.

      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
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
      if ((.NOT.NOTA) && (.NOT.CONJA) && (.NOT.LSAME(TRANSA,'T'))) {
          INFO = 1
      } else if ((.NOT.NOTB) && (.NOT.CONJB) && (.NOT.LSAME(TRANSB,'T'))) {
          INFO = 2
      } else if (M.LT.0) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (K.LT.0) {
          INFO = 5
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 8
      } else if (LDB.LT.MAX(1,NROWB)) {
          INFO = 10
      } else if (LDC.LT.MAX(1,M)) {
          INFO = 13
      }
      if (INFO != 0) {
          xerbla('ZGEMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M == 0) .OR. (N == 0) .OR. (((ALPHA == ZERO).OR. (K == 0)) && (BETA == ONE))) RETURN

      // And when  alpha == zero.

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
          } else if (CONJA) {

            // Form  C := alpha*A**H*B + beta*C.

              for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 100
                          TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
                      } // 100
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 110
              } // 120
          } else {

            // Form  C := alpha*A**T*B + beta*C

              for (J = 1; J <= N; J++) { // 150
                  for (I = 1; I <= M; I++) { // 140
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 130
                          TEMP = TEMP + A(L,I)*B(L,J)
                      } // 130
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
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
                          C(I,J) = ZERO
                      } // 160
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 170
                          C(I,J) = BETA*C(I,J)
                      } // 170
                  }
                  for (L = 1; L <= K; L++) { // 190
                      TEMP = ALPHA*DCONJG(B(J,L))
                      for (I = 1; I <= M; I++) { // 180
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      } // 180
                  } // 190
              } // 200
          } else {

            // Form  C := alpha*A*B**T + beta*C

              for (J = 1; J <= N; J++) { // 250
                  if (BETA == ZERO) {
                      for (I = 1; I <= M; I++) { // 210
                          C(I,J) = ZERO
                      } // 210
                  } else if (BETA != ONE) {
                      for (I = 1; I <= M; I++) { // 220
                          C(I,J) = BETA*C(I,J)
                      } // 220
                  }
                  for (L = 1; L <= K; L++) { // 240
                      TEMP = ALPHA*B(J,L)
                      for (I = 1; I <= M; I++) { // 230
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      } // 230
                  } // 240
              } // 250
          }
      } else if (CONJA) {
          if (CONJB) {

            // Form  C := alpha*A**H*B**H + beta*C.

              for (J = 1; J <= N; J++) { // 280
                  for (I = 1; I <= M; I++) { // 270
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 260
                          TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
                      } // 260
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 270
              } // 280
          } else {

            // Form  C := alpha*A**H*B**T + beta*C

              for (J = 1; J <= N; J++) { // 310
                  for (I = 1; I <= M; I++) { // 300
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 290
                          TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
                      } // 290
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 300
              } // 310
          }
      } else {
          if (CONJB) {

            // Form  C := alpha*A**T*B**H + beta*C

              for (J = 1; J <= N; J++) { // 340
                  for (I = 1; I <= M; I++) { // 330
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 320
                          TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
                      } // 320
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 330
              } // 340
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              for (J = 1; J <= N; J++) { // 370
                  for (I = 1; I <= M; I++) { // 360
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 350
                          TEMP = TEMP + A(L,I)*B(J,L)
                      } // 350
                      if (BETA == ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
                  } // 360
              } // 370
          }
      }

      RETURN

      // End of ZGEMM

      }
