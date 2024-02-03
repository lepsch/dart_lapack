      SUBROUTINE CHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHE, SLAMCH
      // EXTERNAL LSAME, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVHE_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( AIMAG( AFAC( J, J ) ).NE.ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
      } // 10

      // Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // Call CLAVHE_ROOK to form the product D * U' (or D * L' ).

      clavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call CLAVHE_ROOK again to multiply by U (or L ).

      clavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= J - 1; I++) { // 20
               C( I, J ) = C( I, J ) - A( I, J )
            } // 20
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
         } // 30
      } else {
         for (J = 1; J <= N; J++) { // 50
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
            for (I = J + 1; I <= N; I++) { // 40
               C( I, J ) = C( I, J ) - A( I, J )
            } // 40
         } // 50
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS
      }

      RETURN

      // End of CHET01_ROOK

      }
