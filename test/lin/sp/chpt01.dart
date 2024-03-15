      void chpt01(final int UPLO, final int N, final int A, final int AFAC, final Array<int> IPIV_, final Matrix<double> C_, final int LDC, final Array<double> RWORK_, final int RESID,) {
  final IPIV = IPIV_.dim();
  final C = C_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDC, N;
      double               RESID;
      int                IPIV( * );
      double               RWORK( * );
      Complex            A( * ), AFAC( * ), C( LDC, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J, JC;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHE, CLANHP, SLAMCH;
      // EXTERNAL lsame, CLANHE, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAVHP, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHP( '1', UPLO, N, A, RWORK );

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      JC = 1;
      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            if ( AIMAG( AFAC( JC ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
            JC = JC + J + 1;
         } // 10
      } else {
         for (J = 1; J <= N; J++) { // 20
            if ( AIMAG( AFAC( JC ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
            JC = JC + N - J + 1;
         } // 20
      }

      // Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // Call CLAVHP to form the product D * U' (or D * L' ).

      clavhp(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call CLAVHP again to multiply by U ( or L ).

      clavhp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         JC = 0;
         for (J = 1; J <= N; J++) { // 40
            for (I = 1; I <= J - 1; I++) { // 30
               C[I][J] = C( I, J ) - A( JC+I );
            } // 30
            C[J][J] = C( J, J ) - double( A( JC+J ) );
            JC = JC + J;
         } // 40
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 60
            C[J][J] = C( J, J ) - double( A( JC ) );
            for (I = J + 1; I <= N; I++) { // 50
               C[I][J] = C( I, J ) - A( JC+I-J );
            } // 50
            JC = JC + N - J + 1;
         } // 60
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      }