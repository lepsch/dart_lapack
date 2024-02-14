      void zlaqsb(final int UPLO, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final int S, final int SCOND, final int AMAX, final int EQUED,) {
  final AB = AB_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, UPLO;
      int                KD, LDAB, N;
      double             AMAX, SCOND;
      double             S( * );
      Complex         AB( LDAB, * );
      // ..

      double             ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

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

            // Upper triangle of A is stored in band format.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = max( 1, J-KD ); I <= J; I++) { // 10
                  AB[KD+1+I-J][J] = CJ*S( I )*AB( KD+1+I-J, J );
               } // 10
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               for (I = J; I <= min( N, J+KD ); I++) { // 30
                  AB[1+I-J][J] = CJ*S( I )*AB( 1+I-J, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      }
