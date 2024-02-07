      void dsgt01(ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      double             A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double             ANORM, ULP;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYMM

      RESULT[1] = ZERO;
      if (N <= 0) return;

      ULP = dlamch( 'Epsilon' );

      // Compute product of 1-norms of A and Z.

      ANORM = dlansy( '1', UPLO, N, A, LDA, WORK )* dlange( '1', N, M, Z, LDZ, WORK )       IF( ANORM == ZERO ) ANORM = ONE;

      if ( ITYPE == 1 ) {

         // Norm of AZ - BZD

         dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 10
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 10
         dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK, N );

         RESULT[1] = ( dlange( '1', N, M, WORK, N, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 2 ) {

         // Norm of ABZ - ZD

         dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 20
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 20
         dsymm('Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z, LDZ );

         RESULT[1] = ( dlange( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 3 ) {

         // Norm of BAZ - ZD

         dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 30
            dscal(N, D( I ), Z( 1, I ), 1 );
         } // 30
         dsymm('Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z, LDZ );

         RESULT[1] = ( dlange( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );
      }

      return;
      }
