      SUBROUTINE ZGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      double             RESULT( 4 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ), R( LDA, * ), T( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK ), Z( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      COMPLEX*16         CROGUE
      const              CROGUE = ( -1.0D+10, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, RESID, ULP, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZGGQRF, ZHERK, ZLACPY, ZLASET, ZUNGQR, ZUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      zlacpy('Full', N, M, A, LDA, AF, LDA );
      zlacpy('Full', N, P, B, LDB, BF, LDB );

      ANORM = MAX( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      zggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      zlaset('Full', N, N, CROGUE, CROGUE, Q, LDA );
      zlacpy('Lower', N-1, M, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );
      zungqr(N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      zlaset('Full', P, P, CROGUE, CROGUE, Z, LDB );
      if ( N.LE.P ) {
         IF( N.GT.0 .AND. N.LT.P ) CALL ZLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N.GT.1 ) CALL ZLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB )
      } else {
         IF( P.GT.1 ) CALL ZLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB )
      }
      zungrq(P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      zlaset('Full', N, M, CZERO, CZERO, R, LDA );
      zlacpy('Upper', N, M, AF, LDA, R, LDA );

      // Copy T

      zlaset('Full', N, P, CZERO, CZERO, T, LDB );
      if ( N.LE.P ) {
         zlacpy('Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB );
      } else {
         zlacpy('Full', N-P, P, BF, LDB, T, LDB );
         zlacpy('Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB );
      }

      // Compute R - Q'*A

      zgemm('Conjugate transpose', 'No transpose', N, M, N, -CONE, Q, LDA, A, LDA, CONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', N, M, R, LDA, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute T*Z - Q'*B

      zgemm('No Transpose', 'No transpose', N, P, P, CONE, T, LDB, Z, LDB, CZERO, BWK, LDB )       CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, P, N, -CONE, Q, LDA, B, LDB, CONE, BWK, LDB );

      // Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', N, P, BWK, LDB, RWORK )
      if ( BNORM.GT.ZERO ) {
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - Q'*Q

      zlaset('Full', N, N, CZERO, CONE, R, LDA );
      zherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = ZLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP

      // Compute I - Z'*Z

      zlaset('Full', P, P, CZERO, CONE, T, LDB );
      zherk('Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = ZLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP

      RETURN

      // End of ZGQRTS

      }
