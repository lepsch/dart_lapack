      void dgqrts(N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT ) {

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
      // EXTERNAL DGEMM, DGGQRF, DLACPY, DLASET, DORGQR, DORGRQ, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' );
      UNFL = DLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      dlacpy('Full', N, M, A, LDA, AF, LDA );
      dlacpy('Full', N, P, B, LDB, BF, LDB );

      ANORM = max( DLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = max( DLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      dggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      dlacpy('Lower', N-1, M, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );
      dorgqr(N, N, min( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      dlaset('Full', P, P, ROGUE, ROGUE, Z, LDB );
      if ( N <= P ) {
         if (N > 0 && N < P) CALL DLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB );
         IF( N > 1 ) CALL DLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB );
      } else {
         if (P > 1) CALL DLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB );
      }
      dorgrq(P, P, min( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', N, M, ZERO, ZERO, R, LDA );
      dlacpy('Upper', N, M, AF, LDA, R, LDA );

      // Copy T

      dlaset('Full', N, P, ZERO, ZERO, T, LDB );
      if ( N <= P ) {
         dlacpy('Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB );
      } else {
         dlacpy('Full', N-P, P, BF, LDB, T, LDB );
         dlacpy('Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB );
      }

      // Compute R - Q'*A

      dgemm('Transpose', 'No transpose', N, M, N, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( max(M,N)*norm(A)*ULP ) .

      RESID = DLANGE( '1', N, M, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( max( 1, M, N ) ) ) / ANORM ) / ULP;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute T*Z - Q'*B

      dgemm('No Transpose', 'No transpose', N, P, P, ONE, T, LDB, Z, LDB, ZERO, BWK, LDB )       CALL DGEMM( 'Transpose', 'No transpose', N, P, N, -ONE, Q, LDA, B, LDB, ONE, BWK, LDB );

      // Compute norm( T*Z - Q'*B ) / ( max(P,N)*norm(A)*ULP ) .

      RESID = DLANGE( '1', N, P, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT( 2 ) = ( ( RESID / DBLE( max( 1, P, N ) ) ) / BNORM ) / ULP;
      } else {
         RESULT( 2 ) = ZERO;
      }

      // Compute I - Q'*Q

      dlaset('Full', N, N, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

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

      // End of DGQRTS

      }
