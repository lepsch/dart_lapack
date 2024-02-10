      void claqgb(M, N, KL, KU, final Matrix<double> AB, final int LDAB, R, C, ROWCND, COLCND, AMAX, EQUED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED;
      int                KL, KU, LDAB, M, N;
      double               AMAX, COLCND, ROWCND;
      double               C( * ), R( * );
      Complex            AB( LDAB, * );
      // ..

      double               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J;
      double               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

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
               for (I = max( 1, J-KU ); I <= min( M, J+KL ); I++) { // 10
                  AB[KU+1+I-J][J] = CJ*AB( KU+1+I-J, J );
               } // 10
            } // 20
            EQUED = 'C';
         }
      } else if ( COLCND >= THRESH ) {

         // Row scaling, no column scaling

         for (J = 1; J <= N; J++) { // 40
            for (I = max( 1, J-KU ); I <= min( M, J+KL ); I++) { // 30
               AB[KU+1+I-J][J] = R( I )*AB( KU+1+I-J, J );
            } // 30
         } // 40
         EQUED = 'R';
      } else {

         // Row and column scaling

         for (J = 1; J <= N; J++) { // 60
            CJ = C( J );
            for (I = max( 1, J-KU ); I <= min( M, J+KL ); I++) { // 50
               AB[KU+1+I-J][J] = CJ*R( I )*AB( KU+1+I-J, J );
            } // 50
         } // 60
         EQUED = 'B';
      }

      }
