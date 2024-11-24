// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cgbtrs(final int TRANS, final int N, final int KL, final int KU, final int NRHS, final Matrix<double> AB_, final int LDAB, final Array<int> IPIV_, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final AB = AB_.dim();
  final IPIV = IPIV_.dim();
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      int                IPIV( * );
      Complex            AB( LDAB, * ), B( LDB, * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               LNOTI, NOTRAN;
      int                I, J, KD, L, LM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERU, CLACGV, CSWAP, CTBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 ) {
         INFO = -3;
      } else if ( KU < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDAB < ( 2*KL+KU+1 ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CGBTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      KD = KU + KL + 1;
      LNOTI = KL > 0;

      if ( NOTRAN ) {

         // Solve  A*X = B.

         // Solve L*X = B, overwriting B with X.

         // L is represented as a product of permutations and unit lower
         // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
         // where each transformation L(i) is a rank-one modification of
         // the identity matrix.

         if ( LNOTI ) {
            for (J = 1; J <= N - 1; J++) { // 10
               LM = min( KL, N-J );
               L = IPIV( J );
               if (L != J) cswap( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
               cgeru(LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), LDB, B( J+1, 1 ), LDB );
            } // 10
         }

         for (I = 1; I <= NRHS; I++) { // 20

            // Solve U*X = B, overwriting B with X.

            ctbsv('Upper', 'No transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 20

      } else if ( lsame( TRANS, 'T' ) ) {

         // Solve A**T * X = B.

         for (I = 1; I <= NRHS; I++) { // 30

            // Solve U**T * X = B, overwriting B with X.

            ctbsv('Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 30

         // Solve L**T * X = B, overwriting B with X.

         if ( LNOTI ) {
            for (J = N - 1; J >= 1; J--) { // 40
               LM = min( KL, N-J );
               cgemv('Transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               L = IPIV( J );
               if (L != J) cswap( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
            } // 40
         }

      } else {

         // Solve A**H * X = B.

         for (I = 1; I <= NRHS; I++) { // 50

            // Solve U**H * X = B, overwriting B with X.

            ctbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 50

         // Solve L**H * X = B, overwriting B with X.

         if ( LNOTI ) {
            for (J = N - 1; J >= 1; J--) { // 60
               LM = min( KL, N-J );
               clacgv(NRHS, B( J, 1 ), LDB );
               cgemv('Conjugate transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               clacgv(NRHS, B( J, 1 ), LDB );
               L = IPIV( J );
               if (L != J) cswap( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
            } // 60
         }
      }
      }
