import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgetrs(final String TRANS, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      const              ONE = 1.0 ;
      bool               NOTRAN;

      // Test the input parameters.

      INFO.vakue = 0;
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO.vakue = -1;
      } else if ( N < 0 ) {
         INFO.vakue = -2;
      } else if ( NRHS < 0 ) {
         INFO.vakue = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO.vakue = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO.vakue = -8;
      }
      if ( INFO.vakue != 0 ) {
         xerbla('DGETRS', -INFO.vakue );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( NOTRAN ) {

         // Solve A * X = B.

         // Apply row interchanges to the right hand sides.

         dlaswp(NRHS, B, LDB, 1, N, IPIV, 1 );

         // Solve L*X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A**T * X = B.

         // Solve U**T *X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**T *X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Apply row interchanges to the solution vectors.

         dlaswp(NRHS, B, LDB, 1, N, IPIV, -1 );
      }

      }
