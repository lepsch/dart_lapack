// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void ssymv(final int UPLO, final int N, final int ALPHA, final Matrix<double> A_, final int LDA, final int X, final int INCX, final int BETA, final int Y, final int INCY,) {
  final A = A_.dim();

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double ALPHA,BETA;
      int     INCX,INCY,LDA,N;
      String    UPLO;
      double A(LDA,*),X(*),Y(*);
      // ..

      double ONE,ZERO;
      const     ONE=1.0,ZERO=0.0;
      double TEMP1,TEMP2;
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
      // INTRINSIC MAX
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame(UPLO,'U') && !lsame(UPLO,'L')) {
          INFO = 1;
      } else if (N < 0) {
          INFO = 2;
      } else if (LDA < max(1,N)) {
          INFO = 5;
      } else if (INCX == 0) {
          INFO = 7;
      } else if (INCY == 0) {
          INFO = 10;
      }
      if (INFO != 0) {
          xerbla('SSYMV ',INFO);
          return;
      }

      // Quick return if possible.

      if ((N == 0) || ((ALPHA == ZERO) && (BETA == ONE))) return;

      // Set up the start points in  X  and  Y.

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

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      // First form  y := beta*y.

      if (BETA != ONE) {
          if (INCY == 1) {
              if (BETA == ZERO) {
                  for (I = 1; I <= N; I++) { // 10
                      Y[I] = ZERO;
                  } // 10
              } else {
                  for (I = 1; I <= N; I++) { // 20
                      Y[I] = BETA*Y(I);
                  } // 20
              }
          } else {
              IY = KY;
              if (BETA == ZERO) {
                  for (I = 1; I <= N; I++) { // 30
                      Y[IY] = ZERO;
                      IY = IY + INCY;
                  } // 30
              } else {
                  for (I = 1; I <= N; I++) { // 40
                      Y[IY] = BETA*Y(IY);
                      IY = IY + INCY;
                  } // 40
              }
          }
      }
      if (ALPHA == ZERO) return;
      if (lsame(UPLO,'U')) {

         // Form  y  when A is stored in upper triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 60
                  TEMP1 = ALPHA*X(J);
                  TEMP2 = ZERO;
                  for (I = 1; I <= J - 1; I++) { // 50
                      Y[I] = Y(I) + TEMP1*A(I,J);
                      TEMP2 = TEMP2 + A(I,J)*X(I);
                  } // 50
                  Y[J] = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2;
              } // 60
          } else {
              JX = KX;
              JY = KY;
              for (J = 1; J <= N; J++) { // 80
                  TEMP1 = ALPHA*X(JX);
                  TEMP2 = ZERO;
                  IX = KX;
                  IY = KY;
                  for (I = 1; I <= J - 1; I++) { // 70
                      Y[IY] = Y(IY) + TEMP1*A(I,J);
                      TEMP2 = TEMP2 + A(I,J)*X(IX);
                      IX = IX + INCX;
                      IY = IY + INCY;
                  } // 70
                  Y[JY] = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2;
                  JX = JX + INCX;
                  JY = JY + INCY;
              } // 80
          }
      } else {

         // Form  y  when A is stored in lower triangle.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 100
                  TEMP1 = ALPHA*X(J);
                  TEMP2 = ZERO;
                  Y[J] = Y(J) + TEMP1*A(J,J);
                  for (I = J + 1; I <= N; I++) { // 90
                      Y[I] = Y(I) + TEMP1*A(I,J);
                      TEMP2 = TEMP2 + A(I,J)*X(I);
                  } // 90
                  Y[J] = Y(J) + ALPHA*TEMP2;
              } // 100
          } else {
              JX = KX;
              JY = KY;
              for (J = 1; J <= N; J++) { // 120
                  TEMP1 = ALPHA*X(JX);
                  TEMP2 = ZERO;
                  Y[JY] = Y(JY) + TEMP1*A(J,J);
                  IX = JX;
                  IY = JY;
                  for (I = J + 1; I <= N; I++) { // 110
                      IX = IX + INCX;
                      IY = IY + INCY;
                      Y[IY] = Y(IY) + TEMP1*A(I,J);
                      TEMP2 = TEMP2 + A(I,J)*X(IX);
                  } // 110
                  Y[JY] = Y(JY) + ALPHA*TEMP2;
                  JX = JX + INCX;
                  JY = JY + INCY;
              } // 120
          }
      }

      }
