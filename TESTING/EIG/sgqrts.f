      SUBROUTINE SGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), RESULT( 4 ), RWORK( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               ROGUE
      const              ROGUE = -1.0E+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, ULP, UNFL, RESID
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASET, SORGQR, SORGRQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      CALL SLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', N, P, B, LDB, BF, LDB )

      ANORM = MAX( SLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', N, P, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      CALL SGGQRF( N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )

      // Generate the N-by-N matrix Q

      CALL SLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      CALL SLACPY( 'Lower', N-1, M, AF( 2,1 ), LDA, Q( 2,1 ), LDA )
      CALL SORGQR( N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO )

      // Generate the P-by-P matrix Z

      CALL SLASET( 'Full', P, P, ROGUE, ROGUE, Z, LDB )
      if ( N.LE.P ) {
         IF( N.GT.0 .AND. N.LT.P ) CALL SLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N.GT.1 ) CALL SLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB )
      } else {
         IF( P.GT.1) CALL SLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB )
      }
      CALL SORGRQ( P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO )

      // Copy R

      CALL SLASET( 'Full', N, M, ZERO, ZERO, R, LDA )
      CALL SLACPY( 'Upper', N, M, AF, LDA, R, LDA )

      // Copy T

      CALL SLASET( 'Full', N, P, ZERO, ZERO, T, LDB )
      if ( N.LE.P ) {
         CALL SLACPY( 'Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB )
      } else {
         CALL SLACPY( 'Full', N-P, P, BF, LDB, T, LDB )
         CALL SLACPY( 'Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB )
      }

      // Compute R - Q'*A

      CALL SGEMM( 'Transpose', 'No transpose', N, M, N, -ONE, Q, LDA, A, LDA, ONE, R, LDA )

      // Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', N, M, R, LDA, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX(1,M,N) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute T*Z - Q'*B

      CALL SGEMM( 'No Transpose', 'No transpose', N, P, P, ONE, T, LDB, Z, LDB, ZERO, BWK, LDB )       CALL SGEMM( 'Transpose', 'No transpose', N, P, N, -ONE, Q, LDA, B, LDB, ONE, BWK, LDB )

      // Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', N, P, BWK, LDB, RWORK )
      if ( BNORM.GT.ZERO ) {
         RESULT( 2 ) = ( ( RESID / REAL( MAX(1,P,N ) ) )/BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - Q'*Q

      CALL SLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL SSYRK( 'Upper', 'Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP

      // Compute I - Z'*Z

      CALL SLASET( 'Full', P, P, ZERO, ONE, T, LDB )
      CALL SSYRK( 'Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = SLANSY( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP

      RETURN

      // End of SGQRTS

      }
