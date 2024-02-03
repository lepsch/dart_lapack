      SUBROUTINE ZGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )

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
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      COMPLEX*16         CROGUE
      const              CROGUE = ( -1.0e+10, 0.0 ) ;
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
      // EXTERNAL ZGEMM, ZGGRQF, ZHERK, ZLACPY, ZLASET, ZUNGQR, ZUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      zlacpy('Full', M, N, A, LDA, AF, LDA );
      zlacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', P, N, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      zggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      zlaset('Full', N, N, CROGUE, CROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) CALL ZLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M > 1 ) CALL ZLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) CALL ZLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }
      zungrq(N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      zlaset('Full', P, P, CROGUE, CROGUE, Z, LDB );
      if (P > 1) CALL ZLACPY( 'Lower', P-1, N, BF( 2, 1 ), LDB, Z( 2, 1 ), LDB );
      zungqr(P, P, MIN( P, N ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      zlaset('Full', M, N, CZERO, CZERO, R, LDA );
      if ( M <= N ) {
         zlacpy('Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         zlacpy('Full', M-N, N, AF, LDA, R, LDA );
         zlacpy('Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Copy T

      zlaset('Full', P, N, CZERO, CZERO, T, LDB );
      zlacpy('Upper', P, N, BF, LDB, T, LDB );

      // Compute R - A*Q'

      zgemm('No transpose', 'Conjugate transpose', M, N, N, -CONE, A, LDA, Q, LDA, CONE, R, LDA );

      // Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', M, N, R, LDA, RWORK )
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute T*Q - Z'*B

      zgemm('Conjugate transpose', 'No transpose', P, N, P, CONE, Z, LDB, B, LDB, CZERO, BWK, LDB )       CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, CONE, T, LDB, Q, LDA, -CONE, BWK, LDB );

      // Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', P, N, BWK, LDB, RWORK )
      if ( BNORM > ZERO ) {
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, M ) ) ) / BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - Q*Q'

      zlaset('Full', N, N, CZERO, CONE, R, LDA );
      zherk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

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

      // End of ZGRQTS

      }
