      void zlag2c(M, N, A, LDA, SA, LDSA, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            SA( LDSA, * );
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
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..

      RMAX = SLAMCH( 'O' );
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( ( DBLE( A( I, J ) ) < -RMAX ) || ( DBLE( A( I, J ) ) > RMAX ) || ( DIMAG( A( I, J ) ) < -RMAX ) || ( DIMAG( A( I, J ) ) > RMAX ) ) {
               INFO = 1;
               GO TO 30;
            }
            SA( I, J ) = CMPLX( A( I, J ) );
         } // 10
      } // 20
      INFO = 0;
      } // 30
      return;

      // End of ZLAG2C

      }
