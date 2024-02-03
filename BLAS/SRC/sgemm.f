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
      const     ONE=1.0E+0,ZERO=0.0E+0;
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
      if ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. (.NOT.LSAME(TRANSA,'T'))) {
          INFO = 1
      } else if ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. (.NOT.LSAME(TRANSB,'T'))) {
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
      if (INFO.NE.0) {
          xerbla('SGEMM ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN

      // And if  alpha.eq.zero.

      if (ALPHA.EQ.ZERO) {
          if (BETA.EQ.ZERO) {
              for (J = 1; J <= N; J++) { // 20
                  for (I = 1; I <= M; I++) { // 10
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          } else {
              for (J = 1; J <= N; J++) { // 40
                  for (I = 1; I <= M; I++) { // 30
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          }
          RETURN
      }

      // Start the operations.

      if (NOTB) {
          if (NOTA) {

            // Form  C := alpha*A*B + beta*C.

              for (J = 1; J <= N; J++) { // 90
                  if (BETA.EQ.ZERO) {
                      for (I = 1; I <= M; I++) { // 50
                          C(I,J) = ZERO
   50                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      for (I = 1; I <= M; I++) { // 60
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  }
                  for (L = 1; L <= K; L++) { // 80
                      TEMP = ALPHA*B(L,J)
                      for (I = 1; I <= M; I++) { // 70
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          } else {

            // Form  C := alpha*A**T*B + beta*C

              for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 100
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  110             CONTINUE
  120         CONTINUE
          }
      } else {
          if (NOTA) {

            // Form  C := alpha*A*B**T + beta*C

              for (J = 1; J <= N; J++) { // 170
                  if (BETA.EQ.ZERO) {
                      for (I = 1; I <= M; I++) { // 130
                          C(I,J) = ZERO
  130                 CONTINUE
                  } else if (BETA.NE.ONE) {
                      for (I = 1; I <= M; I++) { // 140
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  }
                  for (L = 1; L <= K; L++) { // 160
                      TEMP = ALPHA*B(J,L)
                      for (I = 1; I <= M; I++) { // 150
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          } else {

            // Form  C := alpha*A**T*B**T + beta*C

              for (J = 1; J <= N; J++) { // 200
                  for (I = 1; I <= M; I++) { // 190
                      TEMP = ZERO
                      for (L = 1; L <= K; L++) { // 180
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      if (BETA.EQ.ZERO) {
                          C(I,J) = ALPHA*TEMP
                      } else {
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      }
  190             CONTINUE
  200         CONTINUE
          }
      }

      RETURN

      // End of SGEMM

      }
