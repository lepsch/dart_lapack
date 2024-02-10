      void cher2(final int UPLO, final int N, final int ALPHA, final int X, final int INCX, final int Y, final int INCY, final int A, final int LDA) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex ALPHA;
      int     INCX,INCY,LDA,N;
      String    UPLO;
      Complex A(LDA,*),X(*),Y(*);
      // ..

      Complex ZERO;
      const     ZERO= (0.0,0.0);
      Complex TEMP1,TEMP2;
      int     I,INFO,IX,IY,J,JX,JY,KX,KY;
      // ..
      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX,REAL
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame(UPLO,'U') && !lsame(UPLO,'L')) {
          INFO = 1;
      } else if (N < 0) {
          INFO = 2;
      } else if (INCX == 0) {
          INFO = 5;
      } else if (INCY == 0) {
          INFO = 7;
      } else if (LDA < max(1,N)) {
          INFO = 9;
      }
      if (INFO != 0) {
          xerbla('CHER2 ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || (ALPHA == ZERO)) return;

      // Set up the start points in X and Y if the increments are not both
      // unity.

      if ((INCX != 1) || (INCY != 1)) {
          if (INCX > 0) {
              KX = 1;
          } else {
              KX = 1 - (N-1)*INCX;
          }
          if (INCY > 0) {
              KY = 1;
          } else {
              KY = 1 - (N-1)*INCY;
          }
          JX = KX;
          JY = KY;
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      if (lsame(UPLO,'U')) {

         // Form  A  when A is stored in the upper triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 20
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*CONJG(Y(J));
                      TEMP2 = CONJG(ALPHA*X(J));
                      for (I = 1; I <= J - 1; I++) { // 10
                          A[I][J] = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2;
                      } // 10
                      A[J][J] = REAL(A(J,J)) + double(X(J)*TEMP1+Y(J)*TEMP2);
                  } else {
                      A[J][J] = double(A(J,J));
                  }
              } // 20
          } else {
              for (J = 1; J <= N; J++) { // 40
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*CONJG(Y(JY));
                      TEMP2 = CONJG(ALPHA*X(JX));
                      IX = KX;
                      IY = KY;
                      for (I = 1; I <= J - 1; I++) { // 30
                          A[I][J] = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2;
                          IX = IX + INCX;
                          IY = IY + INCY;
                      } // 30
                      A[J][J] = REAL(A(J,J)) + double(X(JX)*TEMP1+Y(JY)*TEMP2);
                  } else {
                      A[J][J] = double(A(J,J));
                  }
                  JX = JX + INCX;
                  JY = JY + INCY;
              } // 40
          }
      } else {

         // Form  A  when A is stored in the lower triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 60
                  if ((X(J) != ZERO) || (Y(J) != ZERO)) {
                      TEMP1 = ALPHA*CONJG(Y(J));
                      TEMP2 = CONJG(ALPHA*X(J));
                      A[J][J] = REAL(A(J,J)) + double(X(J)*TEMP1+Y(J)*TEMP2);
                      for (I = J + 1; I <= N; I++) { // 50
                          A[I][J] = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2;
                      } // 50
                  } else {
                      A[J][J] = double(A(J,J));
                  }
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  if ((X(JX) != ZERO) || (Y(JY) != ZERO)) {
                      TEMP1 = ALPHA*CONJG(Y(JY));
                      TEMP2 = CONJG(ALPHA*X(JX));
                      A[J][J] = REAL(A(J,J)) + double(X(JX)*TEMP1+Y(JY)*TEMP2);
                      IX = JX;
                      IY = JY;
                      for (I = J + 1; I <= N; I++) { // 70
                          IX = IX + INCX;
                          IY = IY + INCY;
                          A[I][J] = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2;
                      } // 70
                  } else {
                      A[J][J] = double(A(J,J));
                  }
                  JX = JX + INCX;
                  JY = JY + INCY;
              } // 80
          }
      }

      }
