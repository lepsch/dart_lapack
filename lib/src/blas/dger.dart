      import 'dart:math';

import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dger(final int M,final int N,final double ALPHA,final Array<double> X,final int INCX,final Array<double> Y,final int INCY,final Matrix2d<double> A,final int LDA,) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      // double           ALPHA;
      // int     INCX,INCY,LDA,M,N;
      // ..
      // .. Array Arguments ..
      // double           A[LDA,*],X[*],Y[*];
      // ..

// =====================================================================

      // .. Parameters ..
      // double           ZERO;
      const     ZERO=0.0;
      // ..
      // .. Local Scalars ..
      double           TEMP;
      int     I,INFO,IX,J,JY,KX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
          xerbla('DGER  ',INFO);
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
          for (J = 1; J <= N; J++) {
              if (Y[JY] != ZERO) {
                  TEMP = ALPHA*Y[JY];
                  for (I = 1; I <= M; I++) {
                      A[I,J] = A[I,J] + X[I]*TEMP;
                  }
              }
              JY = JY + INCY;
          }
      } else {
          if (INCX > 0) {
              KX = 1;
          } else {
              KX = 1 - (M-1)*INCX;
          }
          for (J = 1; J <= N; J++) {
              if (Y[JY] != ZERO) {
                  TEMP = ALPHA*Y[JY];
                  IX = KX;
                  for (I = 1; I <= M; I++) {
                      A[I,J] = A[I,J] + X[IX]*TEMP;
                      IX = IX + INCX;
                  }
              }
              JY = JY + INCY;
          }
      }

      return;
      }
