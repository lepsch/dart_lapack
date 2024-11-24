// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double slatm2(final int M, final int N, final int I, final int J, final int KL, final int KU, final int IDIST, final Array<int> ISEED_, final int D, final int IGRADE, final int DL, final int DR, final int IPVTNG, final Array<int> IWORK_, final int SPARSE,) {
  final ISEED = ISEED_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      double               SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      double               D( * ), DL( * ), DR( * );
      // ..


      double               ZERO;
      const              ZERO = 0.0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      double               TEMP;
      // ..

      // .. External Functions ..

      double               SLARAN, SLARND;
      // EXTERNAL SLARAN, SLARND
      // ..

// -----------------------------------------------------------------------



      // Check for I and J in range

      if ( I < 1 || I > M || J < 1 || J > N ) {
         SLATM2 = ZERO;
         return;
      }

      // Check for banding

      if ( J > I+KU || J < I-KL ) {
         SLATM2 = ZERO;
         return;
      }

      // Check for sparsity

      if ( SPARSE > ZERO ) {
         if ( SLARAN( ISEED ) < SPARSE ) {
            SLATM2 = ZERO;
            return;
         }
      }

      // Compute subscripts depending on IPVTNG

      if ( IPVTNG == 0 ) {
         ISUB = I;
         JSUB = J;
      } else if ( IPVTNG == 1 ) {
         ISUB = IWORK( I );
         JSUB = J;
      } else if ( IPVTNG == 2 ) {
         ISUB = I;
         JSUB = IWORK( J );
      } else if ( IPVTNG == 3 ) {
         ISUB = IWORK( I );
         JSUB = IWORK( J );
      }

      // Compute entry and grade it according to IGRADE

      if ( ISUB == JSUB ) {
         TEMP = D( ISUB );
      } else {
         TEMP = SLARND( IDIST, ISEED );
      }
      if ( IGRADE == 1 ) {
         TEMP = TEMP*DL( ISUB );
      } else if ( IGRADE == 2 ) {
         TEMP = TEMP*DR( JSUB );
      } else if ( IGRADE == 3 ) {
         TEMP = TEMP*DL( ISUB )*DR( JSUB );
      } else if ( IGRADE == 4 && ISUB != JSUB ) {
         TEMP = TEMP*DL( ISUB ) / DL( JSUB );
      } else if ( IGRADE == 5 ) {
         TEMP = TEMP*DL( ISUB )*DL( JSUB );
      }
      SLATM2 = TEMP;
      }
