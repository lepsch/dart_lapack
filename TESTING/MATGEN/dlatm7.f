      SUBROUTINE DLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, N, RANK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             COND;
      int                IDIST, INFO, IRSIGN, MODE, N, RANK;
      // ..
      // .. Array Arguments ..
      double             D( * );
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D0 ;
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             HALF;
      const              HALF = 0.5D0 ;
      // ..
      // .. Local Scalars ..
      double             ALPHA, TEMP;
      int                I;
      // ..
      // .. External Functions ..
      double             DLARAN;
      // EXTERNAL DLARAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, EXP, LOG
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters. Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

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
         CALL XERBLA( 'DLATM7', -INFO )
         RETURN
      }

      // Compute D according to COND and MODE

      if ( MODE.NE.0 ) {
         GO TO ( 100, 130, 160, 190, 210, 230 )ABS( MODE )

         // One large D value:

  100    CONTINUE
         DO 110 I = 2, RANK
            D( I ) = ONE / COND
  110    CONTINUE
         DO 120 I = RANK + 1, N
            D( I ) = ZERO
  120    CONTINUE
         D( 1 ) = ONE
         GO TO 240

         // One small D value:

  130    CONTINUE
         DO 140 I = 1, RANK - 1
            D( I ) = ONE
  140    CONTINUE
         DO 150 I = RANK + 1, N
            D( I ) = ZERO
  150    CONTINUE
         D( RANK ) = ONE / COND
         GO TO 240

         // Exponentially distributed D values:

  160    CONTINUE
         D( 1 ) = ONE
         if ( N.GT.1 .AND. RANK.GT.1 ) {
            ALPHA = COND**( -ONE / DBLE( RANK-1 ) )
            DO 170 I = 2, RANK
               D( I ) = ALPHA**( I-1 )
  170       CONTINUE
            DO 180 I = RANK + 1, N
               D( I ) = ZERO
  180       CONTINUE
         }
         GO TO 240

         // Arithmetically distributed D values:

  190    CONTINUE
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / DBLE( N-1 )
            DO 200 I = 2, N
               D( I ) = DBLE( N-I )*ALPHA + TEMP
  200       CONTINUE
         }
         GO TO 240

         // Randomly distributed D values on ( 1/COND , 1):

  210    CONTINUE
         ALPHA = LOG( ONE / COND )
         DO 220 I = 1, N
            D( I ) = EXP( ALPHA*DLARAN( ISEED ) )
  220    CONTINUE
         GO TO 240

         // Randomly distributed D values from IDIST

  230    CONTINUE
         CALL DLARNV( IDIST, ISEED, N, D )

  240    CONTINUE

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.1 ) {
            DO 250 I = 1, N
               TEMP = DLARAN( ISEED )
               IF( TEMP.GT.HALF ) D( I ) = -D( I )
  250       CONTINUE
         }

         // Reverse if MODE < 0

         if ( MODE.LT.0 ) {
            DO 260 I = 1, N / 2
               TEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = TEMP
  260       CONTINUE
         }

      }

      RETURN

      // End of DLATM7

      }
