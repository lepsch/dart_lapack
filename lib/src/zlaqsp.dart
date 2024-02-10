      void zlaqsp(final int UPLO, final int N, final int AP, final int S, final int SCOND, final int AMAX, final int EQUED) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, UPLO;
      int                N;
      double             AMAX, SCOND;
      double             S( * );
      Complex         AP( * );
      // ..

      double             ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J, JC;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH

      // Quick return if possible

      if ( N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = dlamch( 'Safe minimum' ) / dlamch( 'Precision' );
      LARGE = ONE / SMALL;

      if ( SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No equilibration

         EQUED = 'N';
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( lsame( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            JC = 1;
            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = 1; I <= J; I++) { // 10
                  AP[JC+I-1] = CJ*S( I )*AP( JC+I-1 );
               } // 10
               JC = JC + J;
            } // 20
         } else {

            // Lower triangle of A is stored.

            JC = 1;
            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               for (I = J; I <= N; I++) { // 30
                  AP[JC+I-J] = CJ*S( I )*AP( JC+I-J );
               } // 30
               JC = JC + N - J + 1;
            } // 40
         }
         EQUED = 'Y';
      }

      }
