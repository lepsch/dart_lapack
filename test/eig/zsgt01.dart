      void zsgt01(final int ITYPE, final int UPLO, final int N, final int M, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> Z, final int LDZ, final int D, final Array<double> _WORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      double             D( * ), RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I;
      double             ANORM, ULP;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZHEMM

      RESULT[1] = ZERO;
      if (N <= 0) return;

      ULP = dlamch( 'Epsilon' );

      // Compute product of 1-norms of A and Z.

      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )* ZLANGE( '1', N, M, Z, LDZ, RWORK )       IF( ANORM == ZERO ) ANORM = ONE;

      if ( ITYPE == 1 ) {

         // Norm of AZ - BZD

         zhemm('Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 10
            zdscal(N, D( I ), Z( 1, I ), 1 );
         } // 10
         zhemm('Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, -CONE, WORK, N );

         RESULT[1] = ( ZLANGE( '1', N, M, WORK, N, RWORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 2 ) {

         // Norm of ABZ - ZD

         zhemm('Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 20
            zdscal(N, D( I ), Z( 1, I ), 1 );
         } // 20
         zhemm('Left', UPLO, N, M, CONE, A, LDA, WORK, N, -CONE, Z, LDZ );

         RESULT[1] = ( ZLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 3 ) {

         // Norm of BAZ - ZD

         zhemm('Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 30
            zdscal(N, D( I ), Z( 1, I ), 1 );
         } // 30
         zhemm('Left', UPLO, N, M, CONE, B, LDB, WORK, N, -CONE, Z, LDZ );

         RESULT[1] = ( ZLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP );
      }

      }
