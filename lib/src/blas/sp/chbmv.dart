// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void chbmv(final int UPLO, final int N, final int K, final int ALPHA, final Matrix<double> A_, final int LDA, final int X, final int INCX, final int BETA, final int Y, final int INCY,) {
  final A = A_.dim();

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      Complex ALPHA,BETA;
      int     INCX,INCY,K,LDA,N;
      String    UPLO;
      Complex A(LDA,*),X(*),Y(*);
      // ..

      Complex ONE;
      const     ONE= (1.0,0.0);
      Complex ZERO;
      const     ZERO= (0.0,0.0);
      Complex TEMP1,TEMP2;
      int     I,INFO,IX,IY,J,JX,JY,KPLUS1,KX,KY,L;
      // ..
      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG,MAX,MIN,REAL
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame(UPLO,'U') && !lsame(UPLO,'L')) {
          INFO = 1;
      } else if (N < 0) {
          INFO = 2;
      } else if (K < 0) {
          INFO = 3;
      } else if (LDA < (K+1)) {
          INFO = 6;
      } else if (INCX == 0) {
          INFO = 8;
      } else if (INCY == 0) {
          INFO = 11;
      }
      if (INFO != 0) {
          xerbla('CHBMV ',INFO);
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

      // Start the operations. In this version the elements of the array A
      // are accessed sequentially with one pass through A.

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

         // Form  y  when upper triangle of A is stored.

          KPLUS1 = K + 1;
          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 60
                  TEMP1 = ALPHA*X(J);
                  TEMP2 = ZERO;
                  L = KPLUS1 - J;
                  for (I = max(1,J-K); I <= J - 1; I++) { // 50
                      Y[I] = Y(I) + TEMP1*A(L+I,J);
                      TEMP2 = TEMP2 + CONJG(A(L+I,J))*X(I);
                  } // 50
                  Y[J] = Y(J) + TEMP1*double(A(KPLUS1,J)) + ALPHA*TEMP2;
              } // 60
          } else {
              JX = KX;
              JY = KY;
              for (J = 1; J <= N; J++) { // 80
                  TEMP1 = ALPHA*X(JX);
                  TEMP2 = ZERO;
                  IX = KX;
                  IY = KY;
                  L = KPLUS1 - J;
                  for (I = max(1,J-K); I <= J - 1; I++) { // 70
                      Y[IY] = Y(IY) + TEMP1*A(L+I,J);
                      TEMP2 = TEMP2 + CONJG(A(L+I,J))*X(IX);
                      IX = IX + INCX;
                      IY = IY + INCY;
                  } // 70
                  Y[JY] = Y(JY) + TEMP1*double(A(KPLUS1,J)) + ALPHA*TEMP2;
                  JX = JX + INCX;
                  JY = JY + INCY;
                  if (J > K) {
                      KX = KX + INCX;
                      KY = KY + INCY;
                  }
              } // 80
          }
      } else {

         // Form  y  when lower triangle of A is stored.

          if ((INCX == 1) && (INCY == 1)) {
              for (J = 1; J <= N; J++) { // 100
                  TEMP1 = ALPHA*X(J);
                  TEMP2 = ZERO;
                  Y[J] = Y(J) + TEMP1*double(A(1,J));
                  L = 1 - J;
                  for (I = J + 1; I <= min(N,J+K); I++) { // 90
                      Y[I] = Y(I) + TEMP1*A(L+I,J);
                      TEMP2 = TEMP2 + CONJG(A(L+I,J))*X(I);
                  } // 90
                  Y[J] = Y(J) + ALPHA*TEMP2;
              } // 100
          } else {
              JX = KX;
              JY = KY;
              for (J = 1; J <= N; J++) { // 120
                  TEMP1 = ALPHA*X(JX);
                  TEMP2 = ZERO;
                  Y[JY] = Y(JY) + TEMP1*double(A(1,J));
                  L = 1 - J;
                  IX = JX;
                  IY = JY;
                  for (I = J + 1; I <= min(N,J+K); I++) { // 110
                      IX = IX + INCX;
                      IY = IY + INCY;
                      Y[IY] = Y(IY) + TEMP1*A(L+I,J);
                      TEMP2 = TEMP2 + CONJG(A(L+I,J))*X(IX);
                  } // 110
                  Y[JY] = Y(JY) + ALPHA*TEMP2;
                  JX = JX + INCX;
                  JY = JY + INCY;
              } // 120
          }
      }

      }
