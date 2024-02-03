      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY);

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA;
      int     INCX,INCY,LDA,M,N;
      String    TRANS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*);
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16 ONE;
      const     ONE= (1.0,0.0);
      COMPLEX*16 ZERO;
      const     ZERO= (0.0,0.0);
      // ..
      // .. Local Scalars ..
      COMPLEX*16 TEMP;
      int     I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY;
      bool    NOCONJ;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG,MAX
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME(TRANS,'N') && !LSAME(TRANS,'T') && !LSAME(TRANS,'C')) {
          INFO = 1;
      } else if (M < 0) {
          INFO = 2;
      } else if (N < 0) {
          INFO = 3;
      } else if (LDA < max(1,M)) {
          INFO = 6;
      } else if (INCX == 0) {
          INFO = 8;
      } else if (INCY == 0) {
          INFO = 11;
      }
      if (INFO != 0) {
          xerbla('ZGEMV ',INFO);
          return;
      }

      // Quick return if possible.

      if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

      NOCONJ = LSAME(TRANS,'T');

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if (LSAME(TRANS,'N')) {
          LENX = N;
          LENY = M;
      } else {
          LENX = M;
          LENY = N;
      }
      if (INCX > 0) {
          KX = 1;
      } else {
          KX = 1 - (LENX-1)*INCX;
      }
      if (INCY > 0) {
          KY = 1;
      } else {
          KY = 1 - (LENY-1)*INCY;
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      // First form  y := beta*y.

      if (BETA != ONE) {
          if (INCY == 1) {
              if (BETA == ZERO) {
                  for (I = 1; I <= LENY; I++) { // 10
                      Y(I) = ZERO;
                  } // 10
              } else {
                  for (I = 1; I <= LENY; I++) { // 20
                      Y(I) = BETA*Y(I);
                  } // 20
              }
          } else {
              IY = KY;
              if (BETA == ZERO) {
                  for (I = 1; I <= LENY; I++) { // 30
                      Y(IY) = ZERO;
                      IY = IY + INCY;
                  } // 30
              } else {
                  for (I = 1; I <= LENY; I++) { // 40
                      Y(IY) = BETA*Y(IY);
                      IY = IY + INCY;
                  } // 40
              }
          }
      }
      if (ALPHA == ZERO) return;
      if (LSAME(TRANS,'N')) {

         // Form  y := alpha*A*x + y.

          JX = KX;
          if (INCY == 1) {
              for (J = 1; J <= N; J++) { // 60
                  TEMP = ALPHA*X(JX);
                  for (I = 1; I <= M; I++) { // 50
                      Y(I) = Y(I) + TEMP*A(I,J);
                  } // 50
                  JX = JX + INCX;
              } // 60
          } else {
              for (J = 1; J <= N; J++) { // 80
                  TEMP = ALPHA*X(JX);
                  IY = KY;
                  for (I = 1; I <= M; I++) { // 70
                      Y(IY) = Y(IY) + TEMP*A(I,J);
                      IY = IY + INCY;
                  } // 70
                  JX = JX + INCX;
              } // 80
          }
      } else {

         // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.

          JY = KY;
          if (INCX == 1) {
              for (J = 1; J <= N; J++) { // 110
                  TEMP = ZERO;
                  if (NOCONJ) {
                      for (I = 1; I <= M; I++) { // 90
                          TEMP = TEMP + A(I,J)*X(I);
                      } // 90
                  } else {
                      for (I = 1; I <= M; I++) { // 100
                          TEMP = TEMP + DCONJG(A(I,J))*X(I);
                      } // 100
                  }
                  Y(JY) = Y(JY) + ALPHA*TEMP;
                  JY = JY + INCY;
              } // 110
          } else {
              for (J = 1; J <= N; J++) { // 140
                  TEMP = ZERO;
                  IX = KX;
                  if (NOCONJ) {
                      for (I = 1; I <= M; I++) { // 120
                          TEMP = TEMP + A(I,J)*X(IX);
                          IX = IX + INCX;
                      } // 120
                  } else {
                      for (I = 1; I <= M; I++) { // 130
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX);
                          IX = IX + INCX;
                      } // 130
                  }
                  Y(JY) = Y(JY) + ALPHA*TEMP;
                  JY = JY + INCY;
              } // 140
          }
      }

      return;

      // End of ZGEMV

      }
