      void sort03(RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       RC;
      int                INFO, K, LDU, LDV, LWORK, MU, MV, N;
      double               RESULT;
      // ..
      // .. Array Arguments ..
      double               U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IRC, J, LMX;
      double               RES1, RES2, S, ULP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ISAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORT01, XERBLA
      // ..
      // .. Executable Statements ..

      // Check inputs

      INFO = 0;
      if ( lsame( RC, 'R' ) ) {
         IRC = 0;
      } else if ( lsame( RC, 'C' ) ) {
         IRC = 1;
      } else {
         IRC = -1;
      }
      if ( IRC == -1 ) {
         INFO = -1;
      } else if ( MU < 0 ) {
         INFO = -2;
      } else if ( MV < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 || K > max( MU, MV ) ) {
         INFO = -5;
      } else if ( ( IRC == 0 && LDU < max( 1, MU ) ) || ( IRC == 1 && LDU < max( 1, N ) ) ) {
         INFO = -7;
      } else if ( ( IRC == 0 && LDV < max( 1, MV ) ) || ( IRC == 1 && LDV < max( 1, N ) ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SORT03', -INFO );
         return;
      }

      // Initialize result

      RESULT = ZERO;
      if (MU == 0 || MV == 0 || N == 0) return;

      // Machine constants

      ULP = SLAMCH( 'Precision' );

      if ( IRC == 0 ) {

         // Compare rows

         RES1 = ZERO;
         for (I = 1; I <= K; I++) { // 20
            LMX = ISAMAX( N, U( I, 1 ), LDU );
            S = SIGN( ONE, U( I, LMX ) )*SIGN( ONE, V( I, LMX ) );
            for (J = 1; J <= N; J++) { // 10
               RES1 = max( RES1, ABS( U( I, J )-S*V( I, J ) ) );
            } // 10
         } // 20
         RES1 = RES1 / ( REAL( N )*ULP );

         // Compute orthogonality of rows of V.

         sort01('Rows', MV, N, V, LDV, WORK, LWORK, RES2 );

      } else {

         // Compare columns

         RES1 = ZERO;
         for (I = 1; I <= K; I++) { // 40
            LMX = ISAMAX( N, U( 1, I ), 1 );
            S = SIGN( ONE, U( LMX, I ) )*SIGN( ONE, V( LMX, I ) );
            for (J = 1; J <= N; J++) { // 30
               RES1 = max( RES1, ABS( U( J, I )-S*V( J, I ) ) );
            } // 30
         } // 40
         RES1 = RES1 / ( REAL( N )*ULP );

         // Compute orthogonality of columns of V.

         sort01('Columns', N, MV, V, LDV, WORK, LWORK, RES2 );
      }

      RESULT = min( max( RES1, RES2 ), ONE / ULP );
      return;
      }