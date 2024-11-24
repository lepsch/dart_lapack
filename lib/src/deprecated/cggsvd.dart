// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cggsvd(final int JOBU, final int JOBV, final int JOBQ, final int M, final int N, final int P, final int K, final int L, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int ALPHA, final int BETA, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final Matrix<double> Q_, final int LDQ, final Array<double> _WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final A = A_.having();
  final B = B_.having();
  final U = U_.having();
  final V = V_.having();
  final Q = Q_.having();
  final _WORK = _WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      int                IWORK( * );
      double               ALPHA( * ), BETA( * ), RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               WANTQ, WANTU, WANTV;
      int                I, IBND, ISUB, J, NCYCLE;
      double               ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGGSVP, CTGSJA, SCOPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Decode and test the input parameters

      WANTU = lsame( JOBU, 'U' );
      WANTV = lsame( JOBV, 'V' );
      WANTQ = lsame( JOBQ, 'Q' );

      INFO = 0;
      if ( !( WANTU || lsame( JOBU, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( WANTV || lsame( JOBV, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( WANTQ || lsame( JOBQ, 'N' ) ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( P < 0 ) {
         INFO = -6;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -10;
      } else if ( LDB < max( 1, P ) ) {
         INFO = -12;
      } else if ( LDU < 1 || ( WANTU && LDU < M ) ) {
         INFO = -16;
      } else if ( LDV < 1 || ( WANTV && LDV < P ) ) {
         INFO = -18;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -20;
      }
      if ( INFO != 0 ) {
         xerbla('CGGSVD', -INFO );
         return;
      }

      // Compute the Frobenius norm of matrices A and B

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK );
      BNORM = CLANGE( '1', P, N, B, LDB, RWORK );

      // Get machine precision and set up threshold for determining
      // the effective numerical rank of the matrices A and B.

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe Minimum' );
      TOLA = max( M, N )*max( ANORM, UNFL )*ULP;
      TOLB = max( P, N )*max( BNORM, UNFL )*ULP;

      cggsvp(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK( N+1 ), INFO );

      // Compute the GSVD of two upper "triangular" matrices

      ctgsja(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO );

      // Sort the singular values and store the pivot indices in IWORK
      // Copy ALPHA to RWORK, then sort ALPHA in RWORK

      scopy(N, ALPHA, 1, RWORK, 1 );
      IBND = min( L, M-K );
      for (I = 1; I <= IBND; I++) {

         // Scan for largest ALPHA(K+I)

         ISUB = I;
         SMAX = RWORK( K+I );
         for (J = I + 1; J <= IBND; J++) {
            TEMP = RWORK( K+J );
            if ( TEMP > SMAX ) {
               ISUB = J;
               SMAX = TEMP;
            }
         }
         if ( ISUB != I ) {
            RWORK[K+ISUB] = RWORK( K+I );
            RWORK[K+I] = SMAX;
            IWORK[K+I] = K + ISUB;
         } else {
            IWORK[K+I] = K + I;
         }
      }

      }
