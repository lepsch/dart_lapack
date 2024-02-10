      void zlat2c(UPLO, N, final Matrix<double> A, final int LDA, final Matrix<double> SA, final int LDSA, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDSA, N;
      Complex            SA( LDSA, * );
      Complex         A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      bool               UPPER;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, CMPLX
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- bool               lsame;
      // EXTERNAL SLAMCH, lsame

      RMAX = SLAMCH( 'O' );
      UPPER = lsame( UPLO, 'U' );
      if ( UPPER ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( ( (A( I, J )).toDouble() < -RMAX ) || ( (A( I, J )).toDouble() > RMAX ) || ( DIMAG( A( I, J ) ) < -RMAX ) || ( DIMAG( A( I, J ) ) > RMAX ) ) {
                  INFO = 1;
                  GO TO 50;
               }
               SA[I][J] = CMPLX( A( I, J ) );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               if ( ( (A( I, J )).toDouble() < -RMAX ) || ( (A( I, J )).toDouble() > RMAX ) || ( DIMAG( A( I, J ) ) < -RMAX ) || ( DIMAG( A( I, J ) ) > RMAX ) ) {
                  INFO = 1;
                  GO TO 50;
               }
               SA[I][J] = CMPLX( A( I, J ) );
            } // 30
         } // 40
      }
      } // 50

      }
