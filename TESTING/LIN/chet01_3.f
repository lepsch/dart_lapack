      SUBROUTINE CHET01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID )

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
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * )
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
      // EXTERNAL CLASET, CLAVHE_ROOK, CSYCONVF_ROOK
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

      // a) Revert to multipliers of L

      csyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO );

      // 1) Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) {
         if ( AIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
      }

      // 2) Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // 3) Call CLAVHE_ROOK to form the product D * U' (or D * L' ).

      clavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 4) Call ZLAVHE_RK again to multiply by U (or L ).

      clavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 5) Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J - 1; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            }
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
         }
      } else {
         for (J = 1; J <= N; J++) {
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
            for (I = J + 1; I <= N; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            }
         }
      }

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS
      }

      // b) Convert to factor of L (or U)

      csyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO );

      RETURN

      // End of CHET01_3

      }
