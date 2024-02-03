      SUBROUTINE DGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ), R( LDA, * ), RESULT( 4 ), RWORK( * ), T( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK ), Z( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, RESID, ULP, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DGGRQF, DLACPY, DLASET, DORGQR, DORGRQ, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' );
      UNFL = DLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      dlacpy('Full', M, N, A, LDA, AF, LDA );
      dlacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = max( DLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = max( DLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      dggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) CALL DLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) CALL DLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) CALL DLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }
      dorgrq(N, N, min( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      dlaset('Full', P, P, ROGUE, ROGUE, Z, LDB );
      if (P > 1) CALL DLACPY( 'Lower', P-1, N, BF( 2, 1 ), LDB, Z( 2, 1 ), LDB );
      dorgqr(P, P, min( P, N ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, R, LDA );
      if ( M <= N ) {
         dlacpy('Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         dlacpy('Full', M-N, N, AF, LDA, R, LDA );
         dlacpy('Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Copy T

      dlaset('Full', P, N, ZERO, ZERO, T, LDB );
      dlacpy('Upper', P, N, BF, LDB, T, LDB );

      // Compute R - A*Q'

      dgemm('No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA );

      // Compute norm( R - A*Q' ) / ( max(M,N)*norm(A)*ULP ) .

      RESID = DLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( max( 1, M, N ) ) ) / ANORM ) / ULP;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute T*Q - Z'*B

      dgemm('Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK, LDB )       CALL DGEMM( 'No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK, LDB );

      // Compute norm( T*Q - Z'*B ) / ( max(P,N)*norm(A)*ULP ) .

      RESID = DLANGE( '1', P, N, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT( 2 ) = ( ( RESID / DBLE( max( 1, P, M ) ) ) / BNORM ) / ULP;
      } else {
         RESULT( 2 ) = ZERO;
      }

      // Compute I - Q*Q'

      dlaset('Full', N, N, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK );
      RESULT( 3 ) = ( RESID / DBLE( max( 1, N ) ) ) / ULP;

      // Compute I - Z'*Z

      dlaset('Full', P, P, ZERO, ONE, T, LDB );
      dsyrk('Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = DLANSY( '1', 'Upper', P, T, LDB, RWORK );
      RESULT( 4 ) = ( RESID / DBLE( max( 1, P ) ) ) / ULP;

      return;

      // End of DGRQTS

      }
