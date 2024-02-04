      void dlag2s(M, N, A, LDA, SA, LDSA, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      double               SA( LDSA, * );
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      RMAX = SLAMCH( 'O' );
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( ( A( I, J ) < -RMAX ) || ( A( I, J ) > RMAX ) ) {
               INFO = 1;
               GO TO 30;
            }
            SA[I, J] = double( A( I, J ) );
         } // 10
      } // 20
      INFO = 0;
      } // 30
      return;
      }
