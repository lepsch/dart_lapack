      double clanhe(final int NORM, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Array<double> WORK_,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM, UPLO;
      int                LDA, N;
      double               WORK( * );
      Complex            A( LDA, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J;
      double               ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, SQRT

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 10
               SUM = ABS( double( A( J, J ) ) );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( double( A( J, J ) ) );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               for (I = J + 1; I <= N; I++) { // 30
                  SUM = ( A( I, J ) ).abs();
                  if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
               } // 30
            } // 40
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO;
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ( A( I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 50
               WORK[J] = SUM + ABS( double( A( J, J ) ) );
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I );
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK[I] = ZERO;
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( double( A( J, J ) ) );
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ( A( I, J ) ).abs();
                  SUM = SUM + ABSA;
                  WORK[I] = WORK( I ) + ABSA;
               } // 90
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 100
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( lsame( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               classq(J-1, A( 1, J ), 1, SCALE, SUM );
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               classq(N-J, A( J+1, J ), 1, SCALE, SUM );
            } // 120
         }
         SUM = 2*SUM;
         for (I = 1; I <= N; I++) { // 130
            if ( double( A( I, I ) ) != ZERO ) {
               ABSA = ABS( double( A( I, I ) ) );
               if ( SCALE < ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2;
                  SCALE = ABSA;
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2;
               }
            }
         } // 130
         VALUE = SCALE*sqrt( SUM );
      }

      CLANHE = VALUE;
      }
