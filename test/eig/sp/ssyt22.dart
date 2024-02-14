      void ssyt22(final int ITYPE, final int UPLO, final int N, final int M, final int KBAND, final Matrix<double> A_, final int LDA, final int D, final int E, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final int TAU, final Array<double> _WORK_, final int RESULT,) {
  final A = A_.dim();
  final U = U_.dim();
  final V = V_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
      double               A( LDA, * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      double               ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANSY;
      // EXTERNAL SLAMCH, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0 || M <= 0) return;

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Precision' );

      // Do Test 1

      // Norm of A:

      ANORM = max( SLANSY( '1', UPLO, N, A, LDA, WORK ), UNFL );

      // Compute error matrix:

      // ITYPE=1: error = U**T A U - S

      ssymm('L', UPLO, N, M, ONE, A, LDA, U, LDU, ZERO, WORK, N );
      NN = N*N;
      NNP1 = NN + 1;
      sgemm('T', 'N', M, M, N, ONE, U, LDU, WORK, N, ZERO, WORK( NNP1 ), N );
      for (J = 1; J <= M; J++) { // 10
         JJ = NN + ( J-1 )*N + J;
         WORK[JJ] = WORK( JJ ) - D( J );
      } // 10
      if ( KBAND == 1 && N > 1 ) {
         for (J = 2; J <= M; J++) { // 20
            JJ1 = NN + ( J-1 )*N + J - 1;
            JJ2 = NN + ( J-2 )*N + J;
            WORK[JJ1] = WORK( JJ1 ) - E( J-1 );
            WORK[JJ2] = WORK( JJ2 ) - E( J-1 );
         } // 20
      }
      WNORM = SLANSY( '1', UPLO, M, WORK( NNP1 ), N, WORK( 1 ) );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( M*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, REAL( M ) ) / ( M*ULP );
         }
      }

      // Do Test 2

      // Compute  U**T U - I

      if (ITYPE == 1) sort01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RESULT( 2 ) );

      }
