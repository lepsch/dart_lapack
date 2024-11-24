// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cgemqr(final int SIDE, final int TRANS, final int M, final int N, final int K, final Matrix<double> A_, final int LDA, final int T, final int TSIZE, final Matrix<double> C_, final int LDC, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final C = C_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, TSIZE, LWORK, LDC;
      Complex            A( LDA, * ), T( * ), C( LDC, * ), WORK( * );
      // ..

// =====================================================================

      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMQRT, CLAMTSQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, MOD

      // Test the input arguments

      LQUERY  = ( LWORK == -1 );
      NOTRAN  = lsame( TRANS, 'N' );
      TRAN    = lsame( TRANS, 'C' );
      LEFT    = lsame( SIDE, 'L' );
      RIGHT   = lsame( SIDE, 'R' );

      MB = INT( T( 2 ) );
      NB = INT( T( 3 ) );
      if ( LEFT ) {
        LW = N * NB;
        MN = M;
      } else {
        LW = MB * NB;
        MN = N;
      }

      MINMNK = min( M, N, K );
      if ( MINMNK == 0 ) {
         LWMIN = 1;
      } else {
         LWMIN = max( 1, LW );
      }

      if ( ( MB > K ) && ( MN > K ) ) {
        if ( ((MN - K) % (MB - K)) == 0 ) {
          NBLCKS = ( MN - K ) / ( MB - K );
        } else {
          NBLCKS = ( MN - K ) / ( MB - K ) + 1;
        }
      } else {
        NBLCKS = 1;
      }

      INFO = 0;
      if ( !LEFT && !RIGHT ) {
        INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
        INFO = -2;
      } else if ( M < 0 ) {
        INFO = -3;
      } else if ( N < 0 ) {
        INFO = -4;
      } else if ( K < 0 || K > MN ) {
        INFO = -5;
      } else if ( LDA < max( 1, MN ) ) {
        INFO = -7;
      } else if ( TSIZE < 5 ) {
        INFO = -9;
      } else if ( LDC < max( 1, M ) ) {
        INFO = -11;
      } else if ( ( LWORK < max( 1, LW ) ) && ( !LQUERY ) ) {
        INFO = -13;
      }

      if ( INFO == 0 ) {
        WORK[1] = SROUNDUP_LWORK( LWMIN );
      }

      if ( INFO != 0 ) {
        xerbla('CGEMQR', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMNK == 0 ) {
        return;
      }

      if( ( LEFT && M <= K ) || ( RIGHT && N <= K ) || ( MB <= K ) || ( MB >= max( M, N, K ) ) ) {
        CALL CGEMQRT( SIDE, TRANS, M, N, K, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, INFO );
      } else {
        clamtsqr(SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, LWORK, INFO );
      }

      WORK[1] = SROUNDUP_LWORK( LWMIN );

      }
