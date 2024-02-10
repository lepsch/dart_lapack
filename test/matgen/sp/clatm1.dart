      void clatm1(MODE, COND, IRSIGN, IDIST, final Array<int> ISEED, D, N, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IDIST, INFO, IRSIGN, MODE, N;
      double               COND;
      int                ISEED( 4 );
      Complex            D( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      int                I;
      double               ALPHA, TEMP;
      Complex            CTEMP;
      // ..
      // .. External Functions ..
      //- REAL               SLARAN;
      //- COMPLEX            CLARND;
      // EXTERNAL SLARAN, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARNV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, EXP, LOG, REAL

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
      } else if ( ( MODE == 6 || MODE == -6 ) && ( IDIST < 1 || IDIST > 4 ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         xerbla('CLATM1', -INFO );
         return;
      }

      // Compute D according to COND and MODE

      if ( MODE != 0 ) {
         GO TO ( 10, 30, 50, 70, 90, 110 )( MODE ).abs();

         // One large D value:

         } // 10
         for (I = 1; I <= N; I++) { // 20
            D[I] = ONE / COND;
         } // 20
         D[1] = ONE;
         GO TO 120;

         // One small D value:

         } // 30
         for (I = 1; I <= N; I++) { // 40
            D[I] = ONE;
         } // 40
         D[N] = ONE / COND;
         GO TO 120;

         // Exponentially distributed D values:

         } // 50
         D[1] = ONE;
         if ( N > 1 ) {
            ALPHA = COND**( -ONE / REAL( N-1 ) );
            for (I = 2; I <= N; I++) { // 60
               D[I] = ALPHA**( I-1 );
            } // 60
         }
         GO TO 120;

         // Arithmetically distributed D values:

         } // 70
         D[1] = ONE;
         if ( N > 1 ) {
            TEMP = ONE / COND;
            ALPHA = ( ONE-TEMP ) / REAL( N-1 );
            for (I = 2; I <= N; I++) { // 80
               D[I] = double( N-I )*ALPHA + TEMP;
            } // 80
         }
         GO TO 120;

         // Randomly distributed D values on ( 1/COND , 1):

         } // 90
         ALPHA = LOG( ONE / COND );
         for (I = 1; I <= N; I++) { // 100
            D[I] = EXP( ALPHA*SLARAN( ISEED ) );
         } // 100
         GO TO 120;

         // Randomly distributed D values from IDIST

         } // 110
         clarnv(IDIST, ISEED, N, D );

         } // 120

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && IRSIGN == 1 ) {
            for (I = 1; I <= N; I++) { // 130
               CTEMP = CLARND( 3, ISEED );
               D[I] = D( I )*( CTEMP / ( CTEMP ).abs() );
            } // 130
         }

         // Reverse if MODE < 0

         if ( MODE < 0 ) {
            for (I = 1; I <= N / 2; I++) { // 140
               CTEMP = D( I );
               D[I] = D( N+1-I );
               D[N+1-I] = CTEMP;
            } // 140
         }

      }

      }
