      void chet22(final int ITYPE, final int UPLO, final int N, final int M, final int KBAND, final Matrix<double> A_, final int LDA, final int D, final int E, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final int TAU, final Array<double> _WORK_, final Array<double> RWORK_, final int RESULT,) {
  final A = A_.dim();
  final U = U_.dim();
  final V = V_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
      double               D( * ), E( * ), RESULT( 2 ), RWORK( * );
      Complex            A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      double               ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- REAL               CLANHE, SLAMCH;
      // EXTERNAL CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHEMM
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

      ANORM = max( CLANHE( '1', UPLO, N, A, LDA, RWORK ), UNFL );

      // Compute error matrix:

      // ITYPE=1: error = U**H A U - S

      chemm('L', UPLO, N, M, CONE, A, LDA, U, LDU, CZERO, WORK, N );
      NN = N*N;
      NNP1 = NN + 1;
      cgemm('C', 'N', M, M, N, CONE, U, LDU, WORK, N, CZERO, WORK( NNP1 ), N );
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
      WNORM = CLANHE( '1', UPLO, M, WORK( NNP1 ), N, RWORK );

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

      // Compute  U**H U - I

      if (ITYPE == 1) cunt01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RWORK, RESULT( 2 ) );

      }