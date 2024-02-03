      SUBROUTINE SGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), RESULT( 4 ), RWORK( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      REAL               ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, ULP, UNFL, RESID;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGGRQF, SLACPY, SLASET, SORGQR, SORGRQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );
      slacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = MAX( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = MAX( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      sggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      slaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) CALL SLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) CALL SLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) CALL SLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }
      sorgrq(N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      slaset('Full', P, P, ROGUE, ROGUE, Z, LDB );
      if (P > 1) CALL SLACPY( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB );
      sorgqr(P, P, MIN( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO );

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

      // Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL(MAX(1,M,N) ) ) / ANORM ) / ULP;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute T*Q - Z'*B

      sgemm('Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK, LDB )       CALL SGEMM( 'No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK, LDB );

      // Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', P, N, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1,P,M ) ) )/BNORM ) / ULP;
      } else {
         RESULT( 2 ) = ZERO;
      }

      // Compute I - Q*Q'

      slaset('Full', N, N, ZERO, ONE, R, LDA );
      ssyrk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK );
      RESULT( 3 ) = ( RESID / REAL( MAX( 1,N ) ) ) / ULP;

      // Compute I - Z'*Z

      slaset('Full', P, P, ZERO, ONE, T, LDB );
      ssyrk('Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = SLANSY( '1', 'Upper', P, T, LDB, RWORK );
      RESULT( 4 ) = ( RESID / REAL( MAX( 1,P ) ) ) / ULP;

      return;

      // End of SGRQTS

      }
