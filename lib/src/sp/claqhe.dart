      void claqhe(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, UPLO;
      int                LDA, N;
      double               AMAX, SCOND;
      double               S( * );
      Complex            A( LDA, * );
      // ..

      double               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J;
      double               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

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

         if ( lsame( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = 1; I <= J - 1; I++) { // 10
                  A[I][J] = CJ*S( I )*A( I, J );
               } // 10
               A[J][J] = CJ*CJ*double( A( J, J ) );
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               A[J][J] = CJ*CJ*double( A( J, J ) );
               for (I = J + 1; I <= N; I++) { // 30
                  A[I][J] = CJ*S( I )*A( I, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      return;
      }
