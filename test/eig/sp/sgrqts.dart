      void sgrqts(final int M, final int P, final int N, final int A, final int AF, final int Q, final int R, final int LDA, final int TAUA, final int B, final int BF, final int Z, final int T, final int BWK, final int LDB, final int TAUB, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, P, N;
      double               A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), RESULT( 4 ), RWORK( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ROGUE;
      const              ROGUE = -1.0e+10 ;
      int                INFO;
      double               ANORM, BNORM, ULP, UNFL, RESID;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGGRQF, SLACPY, SLASET, SORGQR, SORGRQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );
      slacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = max( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = max( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      sggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      slaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) slacpy( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) slacpy( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) slacpy( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }
      sorgrq(N, N, min( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      slaset('Full', P, P, ROGUE, ROGUE, Z, LDB );
      if (P > 1) slacpy( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB );
      sorgqr(P, P, min( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, R, LDA );
      if ( M <= N ) {
         slacpy('Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         slacpy('Full', M-N, N, AF, LDA, R, LDA );
         slacpy('Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Copy T

      slaset('Full', P, N, ZERO, ZERO, T, LDB );
      slacpy('Upper', P, N, BF, LDB, T, LDB );

      // Compute R - A*Q'

      sgemm('No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA );

      // Compute norm( R - A*Q' ) / ( max(M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL(max(1,M,N) ) ) / ANORM ) / ULP;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute T*Q - Z'*B

      sgemm('Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK, LDB )       CALL SGEMM( 'No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK, LDB );

      // Compute norm( T*Q - Z'*B ) / ( max(P,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', P, N, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT[2] = ( ( RESID / REAL( max( 1,P,M ) ) )/BNORM ) / ULP;
      } else {
         RESULT[2] = ZERO;
      }

      // Compute I - Q*Q'

      slaset('Full', N, N, ZERO, ONE, R, LDA );
      ssyrk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK );
      RESULT[3] = ( RESID / REAL( max( 1,N ) ) ) / ULP;

      // Compute I - Z'*Z

      slaset('Full', P, P, ZERO, ONE, T, LDB );
      ssyrk('Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = SLANSY( '1', 'Upper', P, T, LDB, RWORK );
      RESULT[4] = ( RESID / REAL( max( 1,P ) ) ) / ULP;

      }
