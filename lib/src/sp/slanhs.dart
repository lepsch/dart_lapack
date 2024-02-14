      double slanhs(final int NORM, final int N, final Matrix<double> A_, final int LDA, final Array<double> WORK_,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                LDA, N;
      double               A( LDA, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J;
      double               SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= min( N, J+1 ); I++) { // 10
               SUM = ( A( I, J ) ).abs();
               if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
            } // 10
         } // 20
      } else if ( ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO;
            for (I = 1; I <= min( N, J+1 ); I++) { // 30
               SUM = SUM + ( A( I, J ) ).abs();
            } // 30
            if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
         } // 40
      } else if ( lsame( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= N; I++) { // 50
            WORK[I] = ZERO;
         } // 50
         for (J = 1; J <= N; J++) { // 70
            for (I = 1; I <= min( N, J+1 ); I++) { // 60
               WORK[I] = WORK( I ) + ( A( I, J ) ).abs();
            } // 60
         } // 70
         VALUE = ZERO;
         for (I = 1; I <= N; I++) { // 80
            SUM = WORK( I );
            if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
         } // 80
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         for (J = 1; J <= N; J++) { // 90
            slassq(min( N, J+1 ), A( 1, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*sqrt( SUM );
      }

      SLANHS = VALUE;
      }
