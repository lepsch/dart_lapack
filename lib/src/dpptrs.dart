import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dpptrs(UPLO, N, NRHS, AP, final Matrix<double> B, final int LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      double             AP( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DPPTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B where A = U**T * U.

         for (I = 1; I <= NRHS; I++) { // 10

            // Solve U**T *X = B, overwriting B with X.

            dtpsv('Upper', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 );

            // Solve U*X = B, overwriting B with X.

            dtpsv('Upper', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 );
         } // 10
      } else {

         // Solve A*X = B where A = L * L**T.

         for (I = 1; I <= NRHS; I++) { // 20

            // Solve L*Y = B, overwriting B with X.

            dtpsv('Lower', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 );

            // Solve L**T *X = Y, overwriting B with X.

            dtpsv('Lower', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 );
         } // 20
      }

      }
