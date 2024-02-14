      import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/strsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/sp/slaswp.dart';
import 'package:lapack/src/xerbla.dart';

void sgetrs(final String TRANS, final int N, final int NRHS,
  final Matrix<double> A_, final int LDA, final Array<int> IPIV_,
  final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      const              ONE = 1.0 ;
      bool               NOTRAN;

      // Test the input parameters.

      INFO.value = 0;
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO.value = -1;
      } else if ( N < 0 ) {
         INFO.value = -2;
      } else if ( NRHS < 0 ) {
         INFO.value = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO.value = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO.value = -8;
      }
      if ( INFO.value != 0 ) {
         xerbla('SGETRS', -INFO.value );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( NOTRAN ) {

         // Solve A * X = B.

         // Apply row interchanges to the right hand sides.

         slaswp(NRHS, B, LDB, 1, N, IPIV, 1 );

         // Solve L*X = B, overwriting B with X.

         strsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         strsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A**T * X = B.

         // Solve U**T *X = B, overwriting B with X.

         strsm('Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**T *X = B, overwriting B with X.

         strsm('Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Apply row interchanges to the solution vectors.

         slaswp(NRHS, B, LDB, 1, N, IPIV, -1 );
      }

      }
