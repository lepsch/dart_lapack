      void cspt01(UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDC, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * );
      COMPLEX            A( * ), AFAC( * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, JC;
      REAL               ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANSP, CLANSY, SLAMCH;
      // EXTERNAL LSAME, CLANSP, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAVSP, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANSP( '1', UPLO, N, A, RWORK );

      // Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // Call CLAVSP to form the product D * U' (or D * L' ).

      clavsp(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call CLAVSP again to multiply by U ( or L ).

      clavsp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         JC = 0;
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C( I, J ) = C( I, J ) - A( JC+I );
            } // 10
            JC = JC + J;
         } // 20
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C( I, J ) = C( I, J ) - A( JC+I-J );
            } // 30
            JC = JC + N - J + 1;
         } // 40
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS;
      }

      return;

      // End of CSPT01

      }
