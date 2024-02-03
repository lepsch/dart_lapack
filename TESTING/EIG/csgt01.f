      SUBROUTINE CSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ANORM, ULP
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEMM, CSSCAL
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      IF( N.LE.0 ) RETURN

      ULP = SLAMCH( 'Epsilon' )

      // Compute product of 1-norms of A and Z.

      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )* CLANGE( '1', N, M, Z, LDZ, RWORK )       IF( ANORM.EQ.ZERO ) ANORM = ONE

      if ( ITYPE.EQ.1 ) {

         // Norm of AZ - BZD

         chemm('Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 10
            csscal(N, D( I ), Z( 1, I ), 1 );
         } // 10
         chemm('Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, -CONE, WORK, N );

         RESULT( 1 ) = ( CLANGE( '1', N, M, WORK, N, RWORK ) / ANORM ) / ( N*ULP )

      } else if ( ITYPE.EQ.2 ) {

         // Norm of ABZ - ZD

         chemm('Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 20
            csscal(N, D( I ), Z( 1, I ), 1 );
         } // 20
         chemm('Left', UPLO, N, M, CONE, A, LDA, WORK, N, -CONE, Z, LDZ );

         RESULT( 1 ) = ( CLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )

      } else if ( ITYPE.EQ.3 ) {

         // Norm of BAZ - ZD

         chemm('Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 30
            csscal(N, D( I ), Z( 1, I ), 1 );
         } // 30
         chemm('Left', UPLO, N, M, CONE, B, LDB, WORK, N, -CONE, Z, LDZ );

         RESULT( 1 ) = ( CLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )
      }

      RETURN

      // End of CSGT01

      }
