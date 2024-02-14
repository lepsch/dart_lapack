      void ssgt01(final int ITYPE, final int UPLO, final int N, final int M, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Matrix<double> Z_, final int LDZ, final int D, final Array<double> _WORK_, final int RESULT,) {
  final A = A_.dim();
  final B = B_.dim();
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      double               A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double               ANORM, ULP;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYMM

      RESULT[1] = ZERO;
      if (N <= 0) return;

      ULP = SLAMCH( 'Epsilon' );

      // Compute product of 1-norms of A and Z.

      ANORM = SLANSY( '1', UPLO, N, A, LDA, WORK )* SLANGE( '1', N, M, Z, LDZ, WORK )       IF( ANORM == ZERO ) ANORM = ONE;

      if ( ITYPE == 1 ) {

         // Norm of AZ - BZD

         ssymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 10
            sscal(N, D( I ), Z( 1, I ), 1 );
         } // 10
         ssymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK, N );

         RESULT[1] = ( SLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 2 ) {

         // Norm of ABZ - ZD

         ssymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 20
            sscal(N, D( I ), Z( 1, I ), 1 );
         } // 20
         ssymm('Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z, LDZ );

         RESULT[1] = ( SLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );

      } else if ( ITYPE == 3 ) {

         // Norm of BAZ - ZD

         ssymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N );
         for (I = 1; I <= M; I++) { // 30
            sscal(N, D( I ), Z( 1, I ), 1 );
         } // 30
         ssymm('Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z, LDZ );

         RESULT[1] = ( SLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP );
      }

      }
