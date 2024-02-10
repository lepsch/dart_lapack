import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dspsv(UPLO, N, NRHS, AP, IPIV, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      int                IPIV( * );
      double             AP( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSPTRF, DSPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DSPSV ', -INFO );
         return;
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      dsptrf(UPLO, N, AP, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dsptrs(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO );

      }
      }
