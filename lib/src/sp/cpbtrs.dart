// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cpbtrs(final int UPLO, final int N, final int KD, final int NRHS, final Matrix<double> AB_, final int LDAB, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final AB = AB_.dim();
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      Complex            AB( LDAB, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTBSV, XERBLA
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
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD+1 ) {
         INFO = -6;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('CPBTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      if ( UPPER ) {

         // Solve A*X = B where A = U**H *U.

         for (J = 1; J <= NRHS; J++) { // 10

            // Solve U**H *X = B, overwriting B with X.

            ctbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve U*X = B, overwriting B with X.

            ctbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 10
      } else {

         // Solve A*X = B where A = L*L**H.

         for (J = 1; J <= NRHS; J++) { // 20

            // Solve L*X = B, overwriting B with X.

            ctbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );

            // Solve L**H *X = B, overwriting B with X.

            ctbsv('Lower', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB, B( 1, J ), 1 );
         } // 20
      }

      }
