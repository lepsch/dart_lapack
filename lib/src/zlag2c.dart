      void zlag2c(M, N, final Matrix<double> A, final int LDA, final Matrix<double> SA, final int LDSA, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDSA, M, N;
      Complex            SA( LDSA, * );
      Complex         A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, CMPLX
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH

      RMAX = SLAMCH( 'O' );
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( ( (A( I, J )).toDouble() < -RMAX ) || ( (A( I, J )).toDouble() > RMAX ) || ( DIMAG( A( I, J ) ) < -RMAX ) || ( DIMAG( A( I, J ) ) > RMAX ) ) {
               INFO = 1;
               GO TO 30;
            }
            SA[I][J] = CMPLX( A( I, J ) );
         } // 10
      } // 20
      INFO = 0;
      } // 30
      }
