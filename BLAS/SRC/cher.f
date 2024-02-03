      SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL ALPHA
      int     INCX,LDA,N;
      String    UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
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
      // INTRINSIC CONJG,MAX,REAL
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') && .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (N < 0) {
          INFO = 2
      } else if (INCX == 0) {
          INFO = 5
      } else if (LDA < MAX(1,N)) {
          INFO = 7
      }
      if (INFO != 0) {
          xerbla('CHER  ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((N == 0) || (ALPHA == REAL(ZERO))) RETURN

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
                      TEMP = ALPHA*CONJG(X(J))
                      for (I = 1; I <= J - 1; I++) { // 10
                          A(I,J) = A(I,J) + X(I)*TEMP
                      } // 10
                      A(J,J) = REAL(A(J,J)) + REAL(X(J)*TEMP)
                  } else {
                      A(J,J) = REAL(A(J,J))
                  }
              } // 20
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 40
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      IX = KX
                      for (I = 1; I <= J - 1; I++) { // 30
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
                      } // 30
                      A(J,J) = REAL(A(J,J)) + REAL(X(JX)*TEMP)
                  } else {
                      A(J,J) = REAL(A(J,J))
                  }
                  JX = JX + INCX
              } // 40
          }
      } else {

         // Form  A  when A is stored in lower triangle.

          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 60
                  if (X(J) != ZERO) {
                      TEMP = ALPHA*CONJG(X(J))
                      A(J,J) = REAL(A(J,J)) + REAL(TEMP*X(J))
                      for (I = J + 1; I <= N; I++) { // 50
                          A(I,J) = A(I,J) + X(I)*TEMP
                      } // 50
                  } else {
                      A(J,J) = REAL(A(J,J))
                  }
              } // 60
          } else {
              JX = KX
              for (J = 1; J <= N; J++) { // 80
                  if (X(JX) != ZERO) {
                      TEMP = ALPHA*CONJG(X(JX))
                      A(J,J) = REAL(A(J,J)) + REAL(TEMP*X(JX))
                      IX = JX
                      for (I = J + 1; I <= N; I++) { // 70
                          IX = IX + INCX
                          A(I,J) = A(I,J) + X(IX)*TEMP
                      } // 70
                  } else {
                      A(J,J) = REAL(A(J,J))
                  }
                  JX = JX + INCX
              } // 80
          }
      }

      RETURN

      // End of CHER

      }
