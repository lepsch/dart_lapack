      void dsyt22(ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      double             ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANSY;
      // EXTERNAL DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DORT01, DSYMM
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

      ANORM = max( DLANSY( '1', UPLO, N, A, LDA, WORK ), UNFL );

      // Compute error matrix:

      // ITYPE=1: error = U**T A U - S

      dsymm('L', UPLO, N, M, ONE, A, LDA, U, LDU, ZERO, WORK, N );
      NN = N*N;
      NNP1 = NN + 1;
      dgemm('T', 'N', M, M, N, ONE, U, LDU, WORK, N, ZERO, WORK( NNP1 ), N );
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
      WNORM = DLANSY( '1', UPLO, M, WORK( NNP1 ), N, WORK( 1 ) );

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

      // Compute  U**T U - I

      if (ITYPE == 1) CALL DORT01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RESULT( 2 ) );

      return;

      // End of DSYT22

      }
