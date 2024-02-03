      void slaqsy(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                LDA, N;
      REAL               AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), S( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      LARGE = ONE / SMALL;

      if ( SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No equilibration

         EQUED = 'N';
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( LSAME( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = 1; I <= J; I++) { // 10
                  A( I, J ) = CJ*S( I )*A( I, J );
               } // 10
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               for (I = J; I <= N; I++) { // 30
                  A( I, J ) = CJ*S( I )*A( I, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      return;

      // End of SLAQSY

      }
