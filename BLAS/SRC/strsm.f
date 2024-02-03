      SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA
      int     LDA,LDB,M,N;
      String    DIAG,SIDE,TRANSA,UPLO;
      // ..
      // .. Array Arguments ..
      REAL A(LDA,*),B(LDB,*)
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
      int     I,INFO,J,K,NROWA;
      bool    LSIDE,NOUNIT,UPPER;
      // ..
      // .. Parameters ..
      REAL ONE,ZERO
      const     ONE=1.0E+0,ZERO=0.0E+0;
      // ..

      // Test the input parameters.

      LSIDE = LSAME(SIDE,'L')
      if (LSIDE) {
          NROWA = M
      } else {
          NROWA = N
      }
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')

      INFO = 0
      if ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) {
          INFO = 1
      } else if ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) {
          INFO = 2
      } else if ((.NOT.LSAME(TRANSA,'N')) .AND. (.NOT.LSAME(TRANSA,'T')) .AND. (.NOT.LSAME(TRANSA,'C'))) {
          INFO = 3
      } else if ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) {
          INFO = 4
      } else if (M.LT.0) {
          INFO = 5
      } else if (N.LT.0) {
          INFO = 6
      } else if (LDA.LT.MAX(1,NROWA)) {
          INFO = 9
      } else if (LDB.LT.MAX(1,M)) {
          INFO = 11
      }
      if (INFO.NE.0) {
          xerbla('STRSM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF (M.EQ.0 .OR. N.EQ.0) RETURN

      // And when  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          for (J = 1; J <= N; J++) { // 20
              for (I = 1; I <= M; I++) { // 10
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      }

      // Start the operations.

      if (LSIDE) {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*inv( A )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 60
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 30
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      }
                      DO 50 K = M,1,-1
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          }
   50                 CONTINUE
   60             CONTINUE
              } else {
                  for (J = 1; J <= N; J++) { // 100
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 70
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      }
                      for (K = 1; K <= M; K++) { // 90
                          if (B(K,J).NE.ZERO) {
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          }
   90                 CONTINUE
  100             CONTINUE
              }
          } else {

            // Form  B := alpha*inv( A**T )*B.

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 130
                      for (I = 1; I <= M; I++) { // 120
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              } else {
                  for (J = 1; J <= N; J++) { // 160
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              }
          }
      } else {
          if (LSAME(TRANSA,'N')) {

            // Form  B := alpha*B*inv( A ).

              if (UPPER) {
                  for (J = 1; J <= N; J++) { // 210
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 170
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      }
                      DO 190 K = 1,J - 1
                          if (A(K,J).NE.ZERO) {
                              for (I = 1; I <= M; I++) { // 180
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          }
  190                 CONTINUE
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          for (I = 1; I <= M; I++) { // 200
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      }
  210             CONTINUE
              } else {
                  DO 260 J = N,1,-1
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 220
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      }
                      DO 240 K = J + 1,N
                          if (A(K,J).NE.ZERO) {
                              for (I = 1; I <= M; I++) { // 230
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          }
  240                 CONTINUE
                      if (NOUNIT) {
                          TEMP = ONE/A(J,J)
                          for (I = 1; I <= M; I++) { // 250
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      }
  260             CONTINUE
              }
          } else {

            // Form  B := alpha*B*inv( A**T ).

              if (UPPER) {
                  DO 310 K = N,1,-1
                      if (NOUNIT) {
                          TEMP = ONE/A(K,K)
                          for (I = 1; I <= M; I++) { // 270
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      }
                      DO 290 J = 1,K - 1
                          if (A(J,K).NE.ZERO) {
                              TEMP = A(J,K)
                              for (I = 1; I <= M; I++) { // 280
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          }
  290                 CONTINUE
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 300
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      }
  310             CONTINUE
              } else {
                  for (K = 1; K <= N; K++) { // 360
                      if (NOUNIT) {
                          TEMP = ONE/A(K,K)
                          for (I = 1; I <= M; I++) { // 320
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      }
                      DO 340 J = K + 1,N
                          if (A(J,K).NE.ZERO) {
                              TEMP = A(J,K)
                              for (I = 1; I <= M; I++) { // 330
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          }
  340                 CONTINUE
                      if (ALPHA.NE.ONE) {
                          for (I = 1; I <= M; I++) { // 350
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      }
  360             CONTINUE
              }
          }
      }

      RETURN

      // End of STRSM

      }
