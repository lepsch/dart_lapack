      void claqhb(final int UPLO, final int N, final int KD, final Matrix<double> AB, final int LDAB, final int S, final int SCOND, final int AMAX, final int EQUED) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, UPLO;
      int                KD, LDAB, N;
      double               AMAX, SCOND;
      double               S( * );
      Complex            AB( LDAB, * );
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
      // INTRINSIC MAX, MIN, REAL

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

            // Upper triangle of A is stored in band format.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = max( 1, J-KD ); I <= J - 1; I++) { // 10
                  AB[KD+1+I-J][J] = CJ*S( I )*AB( KD+1+I-J, J );
               } // 10
               AB[KD+1][J] = CJ*CJ*double( AB( KD+1, J ) );
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               AB[1][J] = CJ*CJ*double( AB( 1, J ) );
               for (I = J + 1; I <= min( N, J+KD ); I++) { // 30
                  AB[1+I-J][J] = CJ*S( I )*AB( 1+I-J, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      }
