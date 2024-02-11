      void claqge(final int M, final int N, final Matrix<double> A, final int LDA, final int R, final int C, final int ROWCND, final int COLCND, final int AMAX, final int EQUED,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED;
      int                LDA, M, N;
      double               AMAX, COLCND, ROWCND;
      double               C( * ), R( * );
      Complex            A( LDA, * );
      // ..

      double               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J;
      double               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      LARGE = ONE / SMALL;

      if ( ROWCND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No row scaling

         if ( COLCND >= THRESH ) {

            // No column scaling

            EQUED = 'N';
         } else {

            // Column scaling

            for (J = 1; J <= N; J++) { // 20
               CJ = C( J );
               for (I = 1; I <= M; I++) { // 10
                  A[I][J] = CJ*A( I, J );
               } // 10
            } // 20
            EQUED = 'C';
         }
      } else if ( COLCND >= THRESH ) {

         // Row scaling, no column scaling

         for (J = 1; J <= N; J++) { // 40
            for (I = 1; I <= M; I++) { // 30
               A[I][J] = R( I )*A( I, J );
            } // 30
         } // 40
         EQUED = 'R';
      } else {

         // Row and column scaling

         for (J = 1; J <= N; J++) { // 60
            CJ = C( J );
            for (I = 1; I <= M; I++) { // 50
               A[I][J] = CJ*R( I )*A( I, J );
            } // 50
         } // 60
         EQUED = 'B';
      }

      }
