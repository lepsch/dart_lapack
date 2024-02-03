      SUBROUTINE CLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST, INFO, IRSIGN, MODE, N;
      REAL               COND
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            D( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ALPHA, TEMP
      COMPLEX            CTEMP
      // ..
      // .. External Functions ..
      REAL               SLARAN
      COMPLEX            CLARND
      // EXTERNAL SLARAN, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARNV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, EXP, LOG, REAL
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters. Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      if (N == 0) RETURN;

      // Set INFO if an error

      if ( MODE.LT.-6 || MODE.GT.6 ) {
         INFO = -1
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && ( IRSIGN != 0 && IRSIGN != 1 ) ) {
         INFO = -2
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && COND.LT.ONE ) {
         INFO = -3
      } else if ( ( MODE == 6 || MODE == -6 ) && ( IDIST.LT.1 || IDIST.GT.4 ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -7
      }

      if ( INFO != 0 ) {
         xerbla('CLATM1', -INFO );
         RETURN
      }

      // Compute D according to COND and MODE

      if ( MODE != 0 ) {
         GO TO ( 10, 30, 50, 70, 90, 110 )ABS( MODE )

         // One large D value:

         } // 10
         for (I = 1; I <= N; I++) { // 20
            D( I ) = ONE / COND
         } // 20
         D( 1 ) = ONE
         GO TO 120

         // One small D value:

         } // 30
         for (I = 1; I <= N; I++) { // 40
            D( I ) = ONE
         } // 40
         D( N ) = ONE / COND
         GO TO 120

         // Exponentially distributed D values:

         } // 50
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            ALPHA = COND**( -ONE / REAL( N-1 ) )
            for (I = 2; I <= N; I++) { // 60
               D( I ) = ALPHA**( I-1 )
            } // 60
         }
         GO TO 120

         // Arithmetically distributed D values:

         } // 70
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / REAL( N-1 )
            for (I = 2; I <= N; I++) { // 80
               D( I ) = REAL( N-I )*ALPHA + TEMP
            } // 80
         }
         GO TO 120

         // Randomly distributed D values on ( 1/COND , 1):

         } // 90
         ALPHA = LOG( ONE / COND )
         for (I = 1; I <= N; I++) { // 100
            D( I ) = EXP( ALPHA*SLARAN( ISEED ) )
         } // 100
         GO TO 120

         // Randomly distributed D values from IDIST

         } // 110
         clarnv(IDIST, ISEED, N, D );

         } // 120

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && IRSIGN == 1 ) {
            for (I = 1; I <= N; I++) { // 130
               CTEMP = CLARND( 3, ISEED )
               D( I ) = D( I )*( CTEMP / ABS( CTEMP ) )
            } // 130
         }

         // Reverse if MODE < 0

         if ( MODE.LT.0 ) {
            for (I = 1; I <= N / 2; I++) { // 140
               CTEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = CTEMP
            } // 140
         }

      }

      RETURN

      // End of CLATM1

      }
