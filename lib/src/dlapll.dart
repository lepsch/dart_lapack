import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlapll(N, X, INCX, Y, INCY, SSMIN ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, INCY, N;
      double             SSMIN;
      double             X( * ), Y( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             A11, A12, A22, C, SSMAX, TAU;
      // ..
      // .. External Functions ..
      //- double             DDOT;
      // EXTERNAL DDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLARFG, DLAS2

      // Quick return if possible

      if ( N <= 1 ) {
         SSMIN = ZERO;
         return;
      }

      // Compute the QR factorization of the N-by-2 matrix ( X Y )

      dlarfg(N, X( 1 ), X( 1+INCX ), INCX, TAU );
      A11 = X( 1 );
      X[1] = ONE;

      C = -TAU*ddot( N, X, INCX, Y, INCY );
      daxpy(N, C, X, INCX, Y, INCY );

      dlarfg(N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU );

      A12 = Y( 1 );
      A22 = Y( 1+INCY );

      // Compute the SVD of 2-by-2 Upper triangular matrix.

      dlas2(A11, A12, A22, SSMIN, SSMAX );

      return;
      }
