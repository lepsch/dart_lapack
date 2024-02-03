      SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA
      REAL BETA
      int     K,LDA,LDB,LDC,N;
      String    TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
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
      // INTRINSIC CONJG,MAX,REAL
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP1,TEMP2
      int     I,INFO,J,L,NROWA;
      bool    UPPER;
      // ..
      // .. Parameters ..
      REAL ONE
      const     ONE=1.0E+0;
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..

      // Test the input parameters.

      if (LSAME(TRANS,'N')) {
          NROWA = N
      } else {
          NROWA = K
      }
      UPPER = LSAME(UPLO,'U')

      INFO = 0
      if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 1
      } else if ((.NOT.LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'C'))) {
          INFO = 2
      } else if (N.LT.0) {
          INFO = 3
      } else if (K.LT.0) {
          INFO = 4
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 7
      } else if (LDB.LT.MAX(1,NROWA)) {
          INFO = 9
      } else if (LDC.LT.MAX(1,N)) {
          INFO = 12
      }
      if (INFO.NE.0) {
          xerbla('CHER2K',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          if (UPPER) {
              if (BETA.EQ.REAL(ZERO)) {
                  for (J = 1; J <= N; J++) { // 20
                      for (I = 1; I <= J; I++) { // 10
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              } else {
                  for (J = 1; J <= N; J++) { // 40
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = BETA*REAL(C(J,J))
   40             CONTINUE
              }
          } else {
              if (BETA.EQ.REAL(ZERO)) {
                  for (J = 1; J <= N; J++) { // 60
                      for (I = J; I <= N; I++) { // 50
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      C(J,J) = BETA*REAL(C(J,J))
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              }
          }
          RETURN
      }

      // Start the operations.

      if (LSAME(TRANS,'N')) {

         // Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
                    // C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 130
                  if (BETA.EQ.REAL(ZERO)) {
                      for (I = 1; I <= J; I++) { // 90
                          C(I,J) = ZERO
   90                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 100 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                      C(J,J) = BETA*REAL(C(J,J))
                  } else {
                      C(J,J) = REAL(C(J,J))
                  }
                  for (L = 1; L <= K; L++) { // 120
                      if ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) {
                          TEMP1 = ALPHA*CONJG(B(J,L))
                          TEMP2 = CONJG(ALPHA*A(J,L))
                          DO 110 I = 1,J - 1
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  110                     CONTINUE
                          C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2)
                      }
  120             CONTINUE
  130         CONTINUE
          } else {
              for (J = 1; J <= N; J++) { // 180
                  if (BETA.EQ.REAL(ZERO)) {
                      for (I = J; I <= N; I++) { // 140
                          C(I,J) = ZERO
  140                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                      C(J,J) = BETA*REAL(C(J,J))
                  } else {
                      C(J,J) = REAL(C(J,J))
                  }
                  for (L = 1; L <= K; L++) { // 170
                      if ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) {
                          TEMP1 = ALPHA*CONJG(B(J,L))
                          TEMP2 = CONJG(ALPHA*A(J,L))
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  160                     CONTINUE
                          C(J,J) = REAL(C(J,J)) + REAL(A(J,L)*TEMP1+B(J,L)*TEMP2)
                      }
  170             CONTINUE
  180         CONTINUE
          }
      } else {

         // Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
                    // C.

          if (UPPER) {
              for (J = 1; J <= N; J++) { // 210
                  for (I = 1; I <= J; I++) { // 200
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      for (L = 1; L <= K; L++) { // 190
                          TEMP1 = TEMP1 + CONJG(A(L,I))*B(L,J)
                          TEMP2 = TEMP2 + CONJG(B(L,I))*A(L,J)
  190                 CONTINUE
                      if (I.EQ.J) {
                          if (BETA.EQ.REAL(ZERO)) {
                              C(J,J) = REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                          } else {
                              C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                          }
                      } else {
                          if (BETA.EQ.REAL(ZERO)) {
                              C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                          } else {
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                          }
                      }
  200             CONTINUE
  210         CONTINUE
          } else {
              for (J = 1; J <= N; J++) { // 240
                  for (I = J; I <= N; I++) { // 230
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      for (L = 1; L <= K; L++) { // 220
                          TEMP1 = TEMP1 + CONJG(A(L,I))*B(L,J)
                          TEMP2 = TEMP2 + CONJG(B(L,I))*A(L,J)
  220                 CONTINUE
                      if (I.EQ.J) {
                          if (BETA.EQ.REAL(ZERO)) {
                              C(J,J) = REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                          } else {
                              C(J,J) = BETA*REAL(C(J,J)) + REAL(ALPHA*TEMP1+ CONJG(ALPHA)*TEMP2)
                          }
                      } else {
                          if (BETA.EQ.REAL(ZERO)) {
                              C(I,J) = ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                          } else {
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + CONJG(ALPHA)*TEMP2
                          }
                      }
  230             CONTINUE
  240         CONTINUE
          }
      }

      RETURN

      // End of CHER2K

      }
