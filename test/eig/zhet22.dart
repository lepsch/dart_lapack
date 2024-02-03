      void zhet22(ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), RESULT( 2 ), RWORK( * );
      Complex         A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      double             ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANHE;
      // EXTERNAL DLAMCH, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHEMM, ZUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO;
      RESULT( 2 ) = ZERO;
      if (N <= 0 || M <= 0) return;

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Precision' );

      // Do Test 1

      // Norm of A:

      ANORM = max( ZLANHE( '1', UPLO, N, A, LDA, RWORK ), UNFL );

      // Compute error matrix:

      // ITYPE=1: error = U**H A U - S

      zhemm('L', UPLO, N, M, CONE, A, LDA, U, LDU, CZERO, WORK, N );
      NN = N*N;
      NNP1 = NN + 1;
      zgemm('C', 'N', M, M, N, CONE, U, LDU, WORK, N, CZERO, WORK( NNP1 ), N );
      for (J = 1; J <= M; J++) { // 10
         JJ = NN + ( J-1 )*N + J;
         WORK( JJ ) = WORK( JJ ) - D( J );
      } // 10
      if ( KBAND == 1 && N > 1 ) {
         for (J = 2; J <= M; J++) { // 20
            JJ1 = NN + ( J-1 )*N + J - 1;
            JJ2 = NN + ( J-2 )*N + J;
            WORK( JJ1 ) = WORK( JJ1 ) - E( J-1 );
            WORK( JJ2 ) = WORK( JJ2 ) - E( J-1 );
         } // 20
      }
      WNORM = ZLANHE( '1', UPLO, M, WORK( NNP1 ), N, RWORK );

      if ( ANORM > WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT( 1 ) = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*ULP );
         } else {
            RESULT( 1 ) = min( WNORM / ANORM, DBLE( M ) ) / ( M*ULP );
         }
      }

      // Do Test 2

      // Compute  U**H U - I

      if (ITYPE == 1) zunt01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RWORK, RESULT( 2 ) );

      return;
      }
