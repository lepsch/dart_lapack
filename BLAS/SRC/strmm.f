      SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB);

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA;
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),B(LDB,*);
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
      REAL TEMP;
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      REAL ONE,ZERO;
      const     ONE=1.0,ZERO=0.0;
      // ..

      // Test the input parameters.

      LSIDE = LSAME(SIDE,'L');
      if (LSIDE) {
          NROWA = M;
      } else {
          NROWA = N;
      }
      NOUNIT = LSAME(DIAG,'N');
      UPPER = LSAME(UPLO,'U');

      INFO = 0;
      if (( !LSIDE) && ( !LSAME(SIDE,'R'))) {
          INFO = 1;
      } else if (( !UPPER) && ( !LSAME(UPLO,'L'))) {
          INFO = 2;
      } else if (( !LSAME(TRANSA,'N')) && ( !LSAME(TRANSA,'T')) && ( !LSAME(TRANSA,'C'))) {
          INFO = 3;
      } else if (( !LSAME(DIAG,'U')) && ( !LSAME(DIAG,'N'))) {
          INFO = 4;
      } else if (M < 0) {
          INFO = 5;
      } else if (N < 0) {
          INFO = 6;
      } else if (LDA < MAX(1,NROWA)) {
          INFO = 9;
      } else if (LDB < MAX(1,M)) {
          INFO = 11;
      }
      if (INFO != 0) {
          xerbla('STRMM ',INFO);
          return;
      }

      // Quick return if possible.

      if (M == 0 || N == 0) RETURN;

      // And when  alpha == zero.

      if (ALPHA == ZERO) {
          for (J = 1; J <= N; J++) { // 20
              for (I = 1; I <= M; I++) { // 10
                  B(I,J) = ZERO;
              } // 10
          } // 20
          return;
      }

      // Start the operations.

      if (LSIDE) {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*A*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 50
                      for (K = 1; K <= M; K++) { // 40
                          if (B(K,J) != ZERO) {
                              TEMP = ALPHA*B(K,J);
                              for (I = 1; I <= K - 1; I++) { // 30
                                  B(I,J) = B(I,J) + TEMP*A(I,K);
                              } // 30
                              if (NOUNIT) TEMP = TEMP*A(K,K);
                              B(K,J) = TEMP;
                          }
                      } // 40
                  } // 50
              } else {
                  for (J = 1; J <= N; J++) { // 80
                      DO 70 K = M,1,-1;
                          if (B(K,J) != ZERO) {
                              TEMP = ALPHA*B(K,J);
                              B(K,J) = TEMP;
                              if (NOUNIT) B(K,J) = B(K,J)*A(K,K);
                              for (I = K + 1; I <= M; I++) { // 60
                                  B(I,J) = B(I,J) + TEMP*A(I,K);
                              } // 60
                          }
                      } // 70
                  } // 80
              }
          } else {

            // Form  B := alpha*A**T*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 110
                      DO 100 I = M,1,-1;
                          TEMP = B(I,J);
                          if (NOUNIT) TEMP = TEMP*A(I,I);
                          for (K = 1; K <= I - 1; K++) { // 90
                              TEMP = TEMP + A(K,I)*B(K,J);
                          } // 90
                          B(I,J) = ALPHA*TEMP;
                      } // 100
                  } // 110
              } else {
                  for (J = 1; J <= N; J++) { // 140
                      for (I = 1; I <= M; I++) { // 130
                          TEMP = B(I,J);
                          if (NOUNIT) TEMP = TEMP*A(I,I);
                          for (K = I + 1; K <= M; K++) { // 120
                              TEMP = TEMP + A(K,I)*B(K,J);
                          } // 120
                          B(I,J) = ALPHA*TEMP;
                      } // 130
                  } // 140
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*A.

              if (UPPER) {
                  DO 180 J = N,1,-1;
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = 1; I <= M; I++) { // 150
                          B(I,J) = TEMP*B(I,J);
                      } // 150
                      for (K = 1; K <= J - 1; K++) { // 170
                          if (A(K,J) != ZERO) {
                              TEMP = ALPHA*A(K,J);
                              for (I = 1; I <= M; I++) { // 160
                                  B(I,J) = B(I,J) + TEMP*B(I,K);
                              } // 160
                          }
                      } // 170
                  } // 180
              } else {
                  for (J = 1; J <= N; J++) { // 220
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(J,J);
                      for (I = 1; I <= M; I++) { // 190
                          B(I,J) = TEMP*B(I,J);
                      } // 190
                      for (K = J + 1; K <= N; K++) { // 210
                          if (A(K,J) != ZERO) {
                              TEMP = ALPHA*A(K,J);
                              for (I = 1; I <= M; I++) { // 200
                                  B(I,J) = B(I,J) + TEMP*B(I,K);
                              } // 200
                          }
                      } // 210
                  } // 220
              }
          } else {

            // Form  B := alpha*B*A**T.

              if (UPPER) {
                  for (K = 1; K <= N; K++) { // 260
                      for (J = 1; J <= K - 1; J++) { // 240
                          if (A(J,K) != ZERO) {
                              TEMP = ALPHA*A(J,K);
                              for (I = 1; I <= M; I++) { // 230
                                  B(I,J) = B(I,J) + TEMP*B(I,K);
                              } // 230
                          }
                      } // 240
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(K,K);
                      if (TEMP != ONE) {
                          for (I = 1; I <= M; I++) { // 250
                              B(I,K) = TEMP*B(I,K);
                          } // 250
                      }
                  } // 260
              } else {
                  DO 300 K = N,1,-1;
                      for (J = K + 1; J <= N; J++) { // 280
                          if (A(J,K) != ZERO) {
                              TEMP = ALPHA*A(J,K);
                              for (I = 1; I <= M; I++) { // 270
                                  B(I,J) = B(I,J) + TEMP*B(I,K);
                              } // 270
                          }
                      } // 280
                      TEMP = ALPHA;
                      if (NOUNIT) TEMP = TEMP*A(K,K);
                      if (TEMP != ONE) {
                          for (I = 1; I <= M; I++) { // 290
                              B(I,K) = TEMP*B(I,K);
                          } // 290
                      }
                  } // 300
              }
          }
      }

      return;

      // End of STRMM

      }
