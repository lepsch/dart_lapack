// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cpttrs(final int UPLO, final int N, final int NRHS, final int D, final int E, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      double               D( * );
      Complex            B( LDB, * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                IUPLO, J, JB, NB;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CPTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments.

      INFO = 0;
      UPPER = ( UPLO == 'U' || UPLO == 'u' );
      if ( !UPPER && !( UPLO == 'L' || UPLO == 'l' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CPTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = max( 1, ilaenv( 1, 'CPTTRS', UPLO, N, NRHS, -1, -1 ) );
      }

      // Decode UPLO

      if ( UPPER ) {
         IUPLO = 1;
      } else {
         IUPLO = 0;
      }

      if ( NB >= NRHS ) {
         cptts2(IUPLO, N, NRHS, D, E, B, LDB );
      } else {
         for (J = 1; NB < 0 ? J >= NRHS : J <= NRHS; J += NB) { // 10
            JB = min( NRHS-J+1, NB );
            cptts2(IUPLO, N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      }
