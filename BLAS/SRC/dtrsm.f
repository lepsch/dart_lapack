      void dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) {

// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      double           A(LDA,*),B(LDB,*);
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
      double           TEMP;
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      double           ONE,ZERO;
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
      } else if (LDA < max(1,NROWA)) {
          INFO = 9;
      } else if (LDB < max(1,M)) {
          INFO = 11;
      }
      if (INFO != 0) {
          xerbla('DTRSM ',INFO);
          return;
      }

      // Quick return if possible.

      if (M == 0 || N == 0) return;

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

            // Form  B := alpha*inv( A )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 60
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 30
                              B(I,J) = ALPHA*B(I,J);
                          } // 30
                      }
                      DO 50 K = M,1,-1;
                          if (B(K,J) != ZERO) {
                              if (NOUNIT) B(K,J) = B(K,J)/A(K,K);
                              for (I = 1; I <= K - 1; I++) { // 40
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K);
                              } // 40
                          }
                      } // 50
                  } // 60
              } else {
                  for (J = 1; J <= N; J++) { // 100
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 70
                              B(I,J) = ALPHA*B(I,J);
                          } // 70
                      }
                      for (K = 1; K <= M; K++) { // 90
                          if (B(K,J) != ZERO) {
                              if (NOUNIT) B(K,J) = B(K,J)/A(K,K);
                              for (I = K + 1; I <= M; I++) { // 80
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K);
                              } // 80
                          }
                      } // 90
                  } // 100
              }
          } else {

            // Form  B := alpha*inv( A**T )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 130
                      for (I = 1; I <= M; I++) { // 120
                          TEMP = ALPHA*B(I,J);
                          for (K = 1; K <= I - 1; K++) { // 110
                              TEMP = TEMP - A(K,I)*B(K,J);
                          } // 110
                          if (NOUNIT) TEMP = TEMP/A(I,I);
                          B(I,J) = TEMP;
                      } // 120
                  } // 130
              } else {
                  for (J = 1; J <= N; J++) { // 160
                      DO 150 I = M,1,-1;
                          TEMP = ALPHA*B(I,J);
                          for (K = I + 1; K <= M; K++) { // 140
                              TEMP = TEMP - A(K,I)*B(K,J);
                          } // 140
                          if (NOUNIT) TEMP = TEMP/A(I,I);
                          B(I,J) = TEMP;
                      } // 150
                  } // 160
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*inv( A ).

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 210
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 170
                              B(I,J) = ALPHA*B(I,J);
                          } // 170
                      }
                      for (K = 1; K <= J - 1; K++) { // 190
                          if (A(K,J) != ZERO) {
                              for (I = 1; I <= M; I++) { // 180
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K);
                              } // 180
                          }
                      } // 190
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J);
                          for (I = 1; I <= M; I++) { // 200
                              B(I,J) = TEMP*B(I,J);
                          } // 200
                      }
                  } // 210
              } else {
                  DO 260 J = N,1,-1;
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 220
                              B(I,J) = ALPHA*B(I,J);
                          } // 220
                      }
                      for (K = J + 1; K <= N; K++) { // 240
                          if (A(K,J) != ZERO) {
                              for (I = 1; I <= M; I++) { // 230
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K);
                              } // 230
                          }
                      } // 240
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J);
                          for (I = 1; I <= M; I++) { // 250
                              B(I,J) = TEMP*B(I,J);
                          } // 250
                      }
                  } // 260
              }
          } else {

            // Form  B := alpha*B*inv( A**T ).

              if (UPPER) {
                  DO 310 K = N,1,-1;
                      if (NOUNIT) {
                          TEMP = ONE/A(K,K);
                          for (I = 1; I <= M; I++) { // 270
                              B(I,K) = TEMP*B(I,K);
                          } // 270
                      }
                      for (J = 1; J <= K - 1; J++) { // 290
                          if (A(J,K) != ZERO) {
                              TEMP = A(J,K);
                              for (I = 1; I <= M; I++) { // 280
                                  B(I,J) = B(I,J) - TEMP*B(I,K);
                              } // 280
                          }
                      } // 290
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 300
                              B(I,K) = ALPHA*B(I,K);
                          } // 300
                      }
                  } // 310
              } else {
                  for (K = 1; K <= N; K++) { // 360
                      if (NOUNIT) {
                          TEMP = ONE/A(K,K);
                          for (I = 1; I <= M; I++) { // 320
                              B(I,K) = TEMP*B(I,K);
                          } // 320
                      }
                      for (J = K + 1; J <= N; J++) { // 340
                          if (A(J,K) != ZERO) {
                              TEMP = A(J,K);
                              for (I = 1; I <= M; I++) { // 330
                                  B(I,J) = B(I,J) - TEMP*B(I,K);
                              } // 330
                          }
                      } // 340
                      if (ALPHA != ONE) {
                          for (I = 1; I <= M; I++) { // 350
                              B(I,K) = ALPHA*B(I,K);
                          } // 350
                      }
                  } // 360
              }
          }
      }

      return;

      // End of DTRSM

      }
