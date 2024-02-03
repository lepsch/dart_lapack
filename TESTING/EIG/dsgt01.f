      SUBROUTINE DSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, ULP;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYMM
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO;
      if (N <= 0) RETURN;

      ULP = DLAMCH( 'Epsilon' );

      // Compute product of 1-norms of A and Z.

      ANORM = DLANSY( '1', UPLO, N, A, LDA, WORK )* DLANGE( '1', N, M, Z, LDZ, WORK )       IF( ANORM == ZERO ) ANORM = ONE;

      if ( ITYPE == 1 ) {

         // Norm of AZ - BZD

         dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 10
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 10
         dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK, N );

         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 2 ) {

         // Norm of ABZ - ZD

         dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 20
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 20
         dsymm('Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z, LDZ );

         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 3 ) {

         // Norm of BAZ - ZD

         dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 30
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 30
         dsymm('Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z, LDZ );

         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );
      }

      return;

      // End of DSGT01

      }
