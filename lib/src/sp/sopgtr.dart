// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sopgtr(final int UPLO, final int N, final int AP, final int TAU, final Matrix<double> Q_, final int LDQ, final Array<double> _WORK_, final Box<int> INFO,) {
  final Q = Q_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDQ, N;
      double               AP( * ), Q( LDQ, * ), TAU( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               UPPER;
      int                I, IINFO, IJ, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORG2L, SORG2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input arguments

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SOPGTR', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Q was determined by a call to SSPTRD with UPLO = 'U'

         // Unpack the vectors which define the elementary reflectors and
         // set the last row and column of Q equal to those of the unit
         // matrix

         IJ = 2;
         for (J = 1; J <= N - 1; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               Q[I][J] = AP( IJ );
               IJ = IJ + 1;
            } // 10
            IJ = IJ + 2;
            Q[N][J] = ZERO;
         } // 20
         for (I = 1; I <= N - 1; I++) { // 30
            Q[I][N] = ZERO;
         } // 30
         Q[N][N] = ONE;

         // Generate Q(1:n-1,1:n-1)

         sorg2l(N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO );

      } else {

         // Q was determined by a call to SSPTRD with UPLO = 'L'.

         // Unpack the vectors which define the elementary reflectors and
         // set the first row and column of Q equal to those of the unit
         // matrix

         Q[1][1] = ONE;
         for (I = 2; I <= N; I++) { // 40
            Q[I][1] = ZERO;
         } // 40
         IJ = 3;
         for (J = 2; J <= N; J++) { // 60
            Q[1][J] = ZERO;
            for (I = J + 1; I <= N; I++) { // 50
               Q[I][J] = AP( IJ );
               IJ = IJ + 1;
            } // 50
            IJ = IJ + 2;
         } // 60
         if ( N > 1 ) {

            // Generate Q(2:n,2:n)

            sorg2r(N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, IINFO );
         }
      }
      }
