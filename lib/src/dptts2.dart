import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dptts2(N, NRHS, D, E, B, LDB ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDB, N, NRHS;
      double             B( LDB, * ), D( * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL

      // Quick return if possible

      if ( N <= 1 ) {
         if (N == 1) dscal( NRHS, 1.0 / D( 1 ), B, LDB );
         return;
      }

      // Solve A * X = B using the factorization A = L*D*L**T,
      // overwriting each right hand side vector with its solution.

      for (J = 1; J <= NRHS; J++) { // 30

            // Solve L * x = b.

         for (I = 2; I <= N; I++) { // 10
            B[I][J] = B( I, J ) - B( I-1, J )*E( I-1 );
         } // 10

            // Solve D * L**T * x = b.

         B[N][J] = B( N, J ) / D( N );
         for (I = N - 1; I >= 1; I--) { // 20
            B[I][J] = B( I, J ) / D( I ) - B( I+1, J )*E( I );
         } // 20
      } // 30

      }
