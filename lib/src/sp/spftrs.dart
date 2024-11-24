// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void spftrs(final int TRANSR, final int UPLO, final int N, final int NRHS, final int A, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANSR, UPLO;
      int                INFO, LDB, N, NRHS;
      double               A( 0: * ), B( LDB, * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               LOWER, NORMALTRANSR;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, STFSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'T' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('SPFTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // start execution: there are two triangular solves

      if ( LOWER ) {
         stfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
         stfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
      } else {
         stfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
         stfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
      }

      }
