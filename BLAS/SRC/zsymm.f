      SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      int     LDA,LDB,LDC,M,N;
      String    SIDE,UPLO;
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
      // INTRINSIC MAX
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP1,TEMP2
      int     I,INFO,J,K,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      COMPLEX*16 ONE
      const     ONE= (1.0D+0,0.0D+0);
      COMPLEX*16 ZERO
      const     ZERO= (0.0D+0,0.0D+0);
      // ..

      // Set NROWA as the number of rows of A.

      if (LSAME(SIDE,'L')) {
          NROWA = M
      } else {
          NROWA = N
      }
      UPPER = LSAME(UPLO,'U')

      // Test the input parameters.

      INFO = 0
      if ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) {
          INFO = 1
      } else if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 2
      } else if (M.LT.0) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 7
      } else if (LDB.LT.MAX(1,M)) {
          INFO = 9
      } else if (LDC.LT.MAX(1,M)) {
          INFO = 12
      }
      if (INFO.NE.0) {
          xerbla('ZSYMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          if (BETA.EQ.ZERO) {
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

      if (LSAME(SIDE,'L')) {

         // Form  C := alpha*A*B + beta*C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 70
                  for (I = 1; I <= M; I++) { // 60
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
                      } // 50
                      if (BETA.EQ.ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                      }
                  } // 60
              } // 70
          } else {
              for (J = 1; J <= N; J++) { // 100
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
                      } // 80
                      if (BETA.EQ.ZERO) {
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      } else {
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2
                      }
                  } // 90
              } // 100
          }
      } else {

         // Form  C := alpha*B*A + beta*C.

          for (J = 1; J <= N; J++) { // 170
              TEMP1 = ALPHA*A(J,J)
              if (BETA.EQ.ZERO) {
                  for (I = 1; I <= M; I++) { // 110
                      C(I,J) = TEMP1*B(I,J)
                  } // 110
              } else {
                  for (I = 1; I <= M; I++) { // 120
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
                  } // 120
              }
              DO 140 K = 1,J - 1
                  if (UPPER) {
                      TEMP1 = ALPHA*A(K,J)
                  } else {
                      TEMP1 = ALPHA*A(J,K)
                  }
                  for (I = 1; I <= M; I++) { // 130
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
                  } // 130
              } // 140
              DO 160 K = J + 1,N
                  if (UPPER) {
                      TEMP1 = ALPHA*A(J,K)
                  } else {
                      TEMP1 = ALPHA*A(K,J)
                  }
                  for (I = 1; I <= M; I++) { // 150
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
                  } // 150
              } // 160
          } // 170
      }

      RETURN

      // End of ZSYMM

      }
