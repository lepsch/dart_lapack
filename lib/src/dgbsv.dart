import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgbsv(final int N, final int KL, final int KU, final int NRHS, final Matrix<double> AB, final int LDAB, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      int                IPIV( * );
      double             AB( LDAB, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. External Subroutines ..
      // EXTERNAL DGBTRF, DGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( KL < 0 ) {
         INFO = -2;
      } else if ( KU < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDAB < 2*KL+KU+1 ) {
         INFO = -6;
      } else if ( LDB < max( N, 1 ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('DGBSV ', -INFO );
         return;
      }

      // Compute the LU factorization of the band matrix A.

      dgbtrf(N, N, KL, KU, AB, LDAB, IPIV, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dgbtrs('No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO );
      }
      }
