// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void stpmlqt(final int SIDE, final int TRANS, final int M, final int N, final int K, final int L, final int MB, final Matrix<double> V_, final int LDV, final Matrix<double> T_, final int LDT, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Array<double> _WORK_, final Box<int> INFO,) {
  final V = V_.dim();
  final T = T_.dim();
  final A = A_.dim();
  final B = B_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT;
      double   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, NB, LB, KF, LDAQ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, STPRFB
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // .. Test the input arguments ..

      INFO   = 0;
      LEFT   = lsame( SIDE,  'L' );
      RIGHT  = lsame( SIDE,  'R' );
      TRAN   = lsame( TRANS, 'T' );
      NOTRAN = lsame( TRANS, 'N' );

      if ( LEFT ) {
         LDAQ = max( 1, K );
      } else if ( RIGHT ) {
         LDAQ = max( 1, M );
      }
      if ( !LEFT && !RIGHT ) {
         INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 ) {
         INFO = -5;
      } else if ( L < 0 || L > K ) {
         INFO = -6;
      } else if ( MB < 1 || (MB > K && K > 0) ) {
         INFO = -7;
      } else if ( LDV < K ) {
         INFO = -9;
      } else if ( LDT < MB ) {
         INFO = -11;
      } else if ( LDA < LDAQ ) {
         INFO = -13;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -15;
      }

      if ( INFO != 0 ) {
         xerbla('STPMLQT', -INFO );
         return;
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) return;

      if ( LEFT && NOTRAN ) {

         for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
            IB = min( MB, K-I+1 );
            NB = min( M-L+I+IB-1, M );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = 0;
            }
            stprfb('L', 'T', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && TRAN ) {

         for (I = 1; MB < 0 ? I >= K : I <= K; I += MB) {
            IB = min( MB, K-I+1 );
            NB = min( N-L+I+IB-1, N );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = NB-N+L-I+1;
            }
            stprfb('R', 'N', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      } else if ( LEFT && TRAN ) {

         KF = ((K-1)/MB)*MB+1;
         for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
            IB = min( MB, K-I+1 );
            NB = min( M-L+I+IB-1, M );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = 0;
            }
            stprfb('L', 'N', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && NOTRAN ) {

         KF = ((K-1)/MB)*MB+1;
         for (I = KF; -MB < 0 ? I >= 1 : I <= 1; I += -MB) {
            IB = min( MB, K-I+1 );
            NB = min( N-L+I+IB-1, N );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = NB-N+L-I+1;
            }
            stprfb('R', 'T', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      }
