      SUBROUTINE SLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, N, RANK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               COND
      int                IDIST, INFO, IRSIGN, MODE, N, RANK;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               HALF
      const              HALF = 0.5E0 ;
      // ..
      // .. Local Scalars ..
      REAL               ALPHA, TEMP
      int                I;
      // ..
      // .. External Functions ..
      REAL               SLARAN
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

      INFO = 0

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Set INFO if an error

      if ( MODE.LT.-6 .OR. MODE.GT.6 ) {
         INFO = -1
      } else if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. ( IRSIGN.NE.0 .AND. IRSIGN.NE.1 ) ) {
         INFO = -2
      } else if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. COND.LT.ONE ) {
         INFO = -3
      } else if ( ( MODE.EQ.6 .OR. MODE.EQ.-6 ) .AND. ( IDIST.LT.1 .OR. IDIST.GT.3 ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('SLATM7', -INFO );
         RETURN
      }

      // Compute D according to COND and MODE

      if ( MODE.NE.0 ) {
         GO TO ( 100, 130, 160, 190, 210, 230 )ABS( MODE )

         // One large D value:

         } // 100
         for (I = 2; I <= RANK; I++) { // 110
            D( I ) = ONE / COND
         } // 110
         for (I = RANK + 1; I <= N; I++) { // 120
            D( I ) = ZERO
         } // 120
         D( 1 ) = ONE
         GO TO 240

         // One small D value:

         } // 130
         for (I = 1; I <= RANK - 1; I++) { // 140
            D( I ) = ONE
         } // 140
         for (I = RANK + 1; I <= N; I++) { // 150
            D( I ) = ZERO
         } // 150
         D( RANK ) = ONE / COND
         GO TO 240

         // Exponentially distributed D values:

         } // 160
         D( 1 ) = ONE
         if ( N.GT.1  .AND. RANK.GT.1 ) {
            ALPHA = COND**( -ONE / REAL( RANK-1 ) )
            for (I = 2; I <= RANK; I++) { // 170
               D( I ) = ALPHA**( I-1 )
            } // 170
            for (I = RANK + 1; I <= N; I++) { // 180
               D( I ) = ZERO
            } // 180
         }
         GO TO 240

         // Arithmetically distributed D values:

         } // 190
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / REAL( N-1 )
            for (I = 2; I <= N; I++) { // 200
               D( I ) = REAL( N-I )*ALPHA + TEMP
            } // 200
         }
         GO TO 240

         // Randomly distributed D values on ( 1/COND , 1):

         } // 210
         ALPHA = LOG( ONE / COND )
         for (I = 1; I <= N; I++) { // 220
            D( I ) = EXP( ALPHA*SLARAN( ISEED ) )
         } // 220
         GO TO 240

         // Randomly distributed D values from IDIST

         } // 230
         slarnv(IDIST, ISEED, N, D );

         } // 240

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.1 ) {
            for (I = 1; I <= N; I++) { // 250
               TEMP = SLARAN( ISEED )
               if (TEMP.GT.HALF) D( I ) = -D( I );
            } // 250
         }

         // Reverse if MODE < 0

         if ( MODE.LT.0 ) {
            for (I = 1; I <= N / 2; I++) { // 260
               TEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = TEMP
            } // 260
         }

      }

      RETURN

      // End of SLATM7

      }
