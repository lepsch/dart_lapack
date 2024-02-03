      SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           ALPHA;
      int     INCX,LDA,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16 ZERO
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP
      int     I,INFO,IX,J,JX,KX;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE,DCONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0
      if ( !LSAME(UPLO,'U') && !LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N < 0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      } else if (LDA < MAX(1,N)) {
          INFO = 7
      }
      if (INFO != 0) {
          xerbla('ZHER  ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) || (ALPHA == DBLE(ZERO))) RETURN

      // Set the start point in X if the increment is not unity.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX != 1) {
          KX = 1
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      if (LSAME(UPLO,'U')) {

         // Form  A  when A is stored in upper triangle.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 20
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      for (I = 1; I <= J - 1; I++) { // 10
                          A(I,J) = A(I,J) + X(I)*TEMP
                      } // 10
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(J)*TEMP)
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
              } // 20
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      for (I = 1; I <= J - 1; I++) { // 30
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
                      } // 30
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(JX)*TEMP)
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
                  JX = JX + INCX
              } // 40
          }
      } else {

         // Form  A  when A is stored in lower triangle.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 60
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(J))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(J))
                      for (I = J + 1; I <= N; I++) { // 50
                          A(I,J) = A(I,J) + X(I)*TEMP
                      } // 50
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
              } // 60
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*DCONJG(X(JX))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(JX))
                      IX = JX
                      for (I = J + 1; I <= N; I++) { // 70
                          IX = IX + INCX
                          A(I,J) = A(I,J) + X(IX)*TEMP
                      } // 70
                  } else {
                      A(J,J) = DBLE(A(J,J))
                  }
                  JX = JX + INCX
              } // 80
          }
      }

      RETURN

      // End of ZHER

      }
