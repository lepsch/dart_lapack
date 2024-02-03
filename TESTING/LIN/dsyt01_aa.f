      SUBROUTINE DSYT01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DLAVSY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )

      // Initialize C to the tridiagonal matrix T.

      dlaset('Full', N, N, ZERO, ZERO, C, LDC );
      dlacpy('F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 );
      if ( N.GT.1 ) {
         if ( LSAME( UPLO, 'U' ) ) {
            dlacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL DLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
         } else {
            dlacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL DLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
         }

         // Call DTRMM to form the product U' * D (or L * D ).

         if ( LSAME( UPLO, 'U' ) ) {
            dtrmm('Left', UPLO, 'Transpose', 'Unit', N-1, N, ONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC );
         } else {
            dtrmm('Left', UPLO, 'No transpose', 'Unit', N-1, N, ONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC );
         }

         // Call DTRMM again to multiply by U (or L ).

         if ( LSAME( UPLO, 'U' ) ) {
            dtrmm('Right', UPLO, 'No transpose', 'Unit', N, N-1, ONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC );
         } else {
            dtrmm('Right', UPLO, 'Transpose', 'Unit', N, N-1, ONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC );
         }
      }

      // Apply symmetric pivots

      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL DSWAP( N, C( J, 1 ), LDC, C( I, 1 ), LDC )
      END DO
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL DSWAP( N, C( 1, J ), 1, C( 1, I ), 1 )
      END DO


      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = DLANSY( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of DSYT01_AA

      }
