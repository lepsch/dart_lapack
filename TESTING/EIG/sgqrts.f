      SUBROUTINE SGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

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
      // EXTERNAL SGEMM, SLACPY, SLASET, SORGQR, SORGRQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      slacpy('Full', N, M, A, LDA, AF, LDA );
      slacpy('Full', N, P, B, LDB, BF, LDB );

      ANORM = MAX( SLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = MAX( SLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      sggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      slaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      slacpy('Lower', N-1, M, AF( 2,1 ), LDA, Q( 2,1 ), LDA );
      sorgqr(N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      slaset('Full', P, P, ROGUE, ROGUE, Z, LDB );
      if ( N <= P ) {
         if (N > 0 && N < P) CALL SLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N > 1 ) CALL SLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB );
      } else {
         if (P > 1) CALL SLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB );
      }
      sorgrq(P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      slaset('Full', N, M, ZERO, ZERO, R, LDA );
      slacpy('Upper', N, M, AF, LDA, R, LDA );

      // Copy T

      slaset('Full', N, P, ZERO, ZERO, T, LDB );
      if ( N <= P ) {
         slacpy('Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB );
      } else {
         slacpy('Full', N-P, P, BF, LDB, T, LDB );
         slacpy('Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB );
      }

      // Compute R - Q'*A

      sgemm('Transpose', 'No transpose', N, M, N, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', N, M, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX(1,M,N) ) ) / ANORM ) / ULP;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute T*Z - Q'*B

      sgemm('No Transpose', 'No transpose', N, P, P, ONE, T, LDB, Z, LDB, ZERO, BWK, LDB )       CALL SGEMM( 'Transpose', 'No transpose', N, P, N, -ONE, Q, LDA, B, LDB, ONE, BWK, LDB );

      // Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', N, P, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT( 2 ) = ( ( RESID / REAL( MAX(1,P,N ) ) )/BNORM ) / ULP;
      } else {
         RESULT( 2 ) = ZERO;
      }

      // Compute I - Q'*Q

      slaset('Full', N, N, ZERO, ONE, R, LDA );
      ssyrk('Upper', 'Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK );
      RESULT( 3 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP;

      // Compute I - Z'*Z

      slaset('Full', P, P, ZERO, ONE, T, LDB );
      ssyrk('Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = SLANSY( '1', 'Upper', P, T, LDB, RWORK );
      RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP;

      return;

      // End of SGQRTS

      }
