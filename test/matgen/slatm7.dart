      void slatm7(MODE, COND, IRSIGN, IDIST, ISEED, D, N, RANK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               COND;
      int                IDIST, INFO, IRSIGN, MODE, N, RANK;
      // ..
      // .. Array Arguments ..
      double               D( * );
      int                ISEED( 4 );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      double               HALF;
      const              HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      double               ALPHA, TEMP;
      int                I;
      // ..
      // .. External Functions ..
      //- REAL               SLARAN;
      // EXTERNAL SLARAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARNV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, EXP, LOG, REAL
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters. Initialize flags & seed.

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;

      // Set INFO if an error

      if ( MODE < -6 || MODE > 6 ) {
         INFO = -1;
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && ( IRSIGN != 0 && IRSIGN != 1 ) ) {
         INFO = -2;
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && COND < ONE ) {
         INFO = -3;
      } else if ( ( MODE == 6 || MODE == -6 ) && ( IDIST < 1 || IDIST > 3 ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         xerbla('SLATM7', -INFO );
         return;
      }

      // Compute D according to COND and MODE

      if ( MODE != 0 ) {
         GO TO ( 100, 130, 160, 190, 210, 230 )( MODE ).abs();

         // One large D value:

         } // 100
         for (I = 2; I <= RANK; I++) { // 110
            D[I] = ONE / COND;
         } // 110
         for (I = RANK + 1; I <= N; I++) { // 120
            D[I] = ZERO;
         } // 120
         D[1] = ONE;
         GO TO 240;

         // One small D value:

         } // 130
         for (I = 1; I <= RANK - 1; I++) { // 140
            D[I] = ONE;
         } // 140
         for (I = RANK + 1; I <= N; I++) { // 150
            D[I] = ZERO;
         } // 150
         D[RANK] = ONE / COND;
         GO TO 240;

         // Exponentially distributed D values:

         } // 160
         D[1] = ONE;
         if ( N > 1 && RANK > 1 ) {
            ALPHA = COND**( -ONE / REAL( RANK-1 ) );
            for (I = 2; I <= RANK; I++) { // 170
               D[I] = ALPHA**( I-1 );
            } // 170
            for (I = RANK + 1; I <= N; I++) { // 180
               D[I] = ZERO;
            } // 180
         }
         GO TO 240;

         // Arithmetically distributed D values:

         } // 190
         D[1] = ONE;
         if ( N > 1 ) {
            TEMP = ONE / COND;
            ALPHA = ( ONE-TEMP ) / REAL( N-1 );
            for (I = 2; I <= N; I++) { // 200
               D[I] = double( N-I )*ALPHA + TEMP;
            } // 200
         }
         GO TO 240;

         // Randomly distributed D values on ( 1/COND , 1):

         } // 210
         ALPHA = LOG( ONE / COND );
         for (I = 1; I <= N; I++) { // 220
            D[I] = EXP( ALPHA*SLARAN( ISEED ) );
         } // 220
         GO TO 240;

         // Randomly distributed D values from IDIST

         } // 230
         slarnv(IDIST, ISEED, N, D );

         } // 240

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && IRSIGN == 1 ) {
            for (I = 1; I <= N; I++) { // 250
               TEMP = SLARAN( ISEED );
               if (TEMP > HALF) D( I ) = -D( I );
            } // 250
         }

         // Reverse if MODE < 0

         if ( MODE < 0 ) {
            for (I = 1; I <= N / 2; I++) { // 260
               TEMP = D( I );
               D[I] = D( N+1-I );
               D[N+1-I] = TEMP;
            } // 260
         }

      }

      return;
      }