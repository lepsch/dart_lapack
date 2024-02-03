      SUBROUTINE ZLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST, INFO, IRSIGN, MODE, N;
      double             COND;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX*16         D( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ALPHA, TEMP;
      COMPLEX*16         CTEMP
      // ..
      // .. External Functions ..
      double             DLARAN;
      COMPLEX*16         ZLARND
      // EXTERNAL DLARAN, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARNV
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
      } else if ( ( MODE.EQ.6 .OR. MODE.EQ.-6 ) .AND. ( IDIST.LT.1 .OR. IDIST.GT.4 ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('ZLATM1', -INFO );
         RETURN
      }

      // Compute D according to COND and MODE

      if ( MODE.NE.0 ) {
         GO TO ( 10, 30, 50, 70, 90, 110 )ABS( MODE )

         // One large D value:

   10    CONTINUE
         DO 20 I = 1, N
            D( I ) = ONE / COND
   20    CONTINUE
         D( 1 ) = ONE
         GO TO 120

         // One small D value:

   30    CONTINUE
         DO 40 I = 1, N
            D( I ) = ONE
   40    CONTINUE
         D( N ) = ONE / COND
         GO TO 120

         // Exponentially distributed D values:

   50    CONTINUE
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            ALPHA = COND**( -ONE / DBLE( N-1 ) )
            DO 60 I = 2, N
               D( I ) = ALPHA**( I-1 )
   60       CONTINUE
         }
         GO TO 120

         // Arithmetically distributed D values:

   70    CONTINUE
         D( 1 ) = ONE
         if ( N.GT.1 ) {
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / DBLE( N-1 )
            DO 80 I = 2, N
               D( I ) = DBLE( N-I )*ALPHA + TEMP
   80       CONTINUE
         }
         GO TO 120

         // Randomly distributed D values on ( 1/COND , 1):

   90    CONTINUE
         ALPHA = LOG( ONE / COND )
         DO 100 I = 1, N
            D( I ) = EXP( ALPHA*DLARAN( ISEED ) )
  100    CONTINUE
         GO TO 120

         // Randomly distributed D values from IDIST

  110    CONTINUE
         zlarnv(IDIST, ISEED, N, D );

  120    CONTINUE

         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D

         if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.1 ) {
            DO 130 I = 1, N
               CTEMP = ZLARND( 3, ISEED )
               D( I ) = D( I )*( CTEMP / ABS( CTEMP ) )
  130       CONTINUE
         }

         // Reverse if MODE < 0

         if ( MODE.LT.0 ) {
            DO 140 I = 1, N / 2
               CTEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = CTEMP
  140       CONTINUE
         }

      }

      RETURN

      // End of ZLATM1

      }
