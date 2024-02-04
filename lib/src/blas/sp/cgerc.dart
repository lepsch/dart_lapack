      void cgerc(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex ALPHA;
      int     INCX,INCY,LDA,M,N;
      // ..
      // .. Array Arguments ..
      Complex A(LDA,*),X(*),Y(*);
      // ..

// =====================================================================

      // .. Parameters ..
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      Complex TEMP;
      int     I,INFO,IX,J,JY,KX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0;
      if (M < 0) {
          INFO = 1;
      } else if (N < 0) {
          INFO = 2;
      } else if (INCX == 0) {
          INFO = 5;
      } else if (INCY == 0) {
          INFO = 7;
      } else if (LDA < max(1,M)) {
          INFO = 9;
      }
      if (INFO != 0) {
          xerbla('CGERC ',INFO);
          return;
      }

      // Quick return if possible.

      if ((M == 0) || (N == 0) || (ALPHA == ZERO)) return;

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (INCY > 0) {
          JY = 1;
      } else {
          JY = 1 - (N-1)*INCY;
      }
      if (INCX == 1) {
          for (J = 1; J <= N; J++) { // 20
              if (Y(JY) != ZERO) {
                  TEMP = ALPHA*CONJG(Y(JY));
                  for (I = 1; I <= M; I++) { // 10
                      A[I,J] = A(I,J) + X(I)*TEMP;
                  } // 10
              }
              JY = JY + INCY;
          } // 20
      } else {
          if (INCX > 0) {
              KX = 1;
          } else {
              KX = 1 - (M-1)*INCX;
          }
          for (J = 1; J <= N; J++) { // 40
              if (Y(JY) != ZERO) {
                  TEMP = ALPHA*CONJG(Y(JY));
                  IX = KX;
                  for (I = 1; I <= M; I++) { // 30
                      A[I,J] = A(I,J) + X(IX)*TEMP;
                      IX = IX + INCX;
                  } // 30
              }
              JY = JY + INCY;
          } // 40
      }

      return;
      }