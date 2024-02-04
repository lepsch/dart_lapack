      void csyt01_3(UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, CLANSY;
      // EXTERNAL lsame, SLAMCH, CLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVSY_ROOK, CSYCONVF_ROOK
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

      // a) Revert to multipliers of L

      csyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO );

      // 1) Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK );

      // 2) Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // 3) Call ZLAVSY_ROOK to form the product D * U' (or D * L' ).

      clavsy_rook(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 4) Call ZLAVSY_ROOK again to multiply by U (or L ).

      clavsy_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 5) Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               C[I, J] = C( I, J ) - A( I, J );
            }
         }
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               C[I, J] = C( I, J ) - A( I, J );
            }
         }
      }

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }


      // b) Convert to factor of L (or U)

      csyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO );

      return;
      }
