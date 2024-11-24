// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cqlt02(final int M, final int N, final int K, final int A, final int AF, final int Q, final int L, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               RESULT( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNGQL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Quick return if possible

      if ( M == 0 || N == 0 || K == 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );

      // Copy the last k columns of the factorization to the array Q

      claset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < M) clacpy( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, N-K+1 ), LDA );
      IF( K > 1 ) clacpy( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, N-K+2 ), LDA );

      // Generate the last n columns of the matrix Q

     srnamc.SRNAMT = 'CUNGQL';
      cungql(M, N, K, Q, LDA, TAU( N-K+1 ), WORK, LWORK, INFO );

      // Copy L(m-n+1:m,n-k+1:n)

      claset('Full', N, K, CMPLX( ZERO ), CMPLX( ZERO ), L( M-N+1, N-K+1 ), LDA )       CALL CLACPY( 'Lower', K, K, AF( M-K+1, N-K+1 ), LDA, L( M-K+1, N-K+1 ), LDA );

      // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

      cgemm('Conjugate transpose', 'No transpose', N, K, M, CMPLX( -ONE ), Q, LDA, A( 1, N-K+1 ), LDA, CMPLX( ONE ), L( M-N+1, N-K+1 ), LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, K, A( 1, N-K+1 ), LDA, RWORK );
      RESID = CLANGE( '1', N, K, L( M-N+1, N-K+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), L, LDA );
      cherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      }
