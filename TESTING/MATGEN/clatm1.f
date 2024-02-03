      SUBROUTINE CLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                IDIST, INFO, IRSIGN, MODE, N;
      REAL               COND
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            D( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
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
*
      // Decode and Test the input parameters. Initialize flags & seed.
*
      INFO = 0
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Set INFO if an error
*
      IF( MODE.LT.-6 .OR. MODE.GT.6 ) THEN
         INFO = -1
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. ( IRSIGN.NE.0 .AND. IRSIGN.NE.1 ) ) THEN
         INFO = -2
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. COND.LT.ONE ) THEN
         INFO = -3
      ELSE IF( ( MODE.EQ.6 .OR. MODE.EQ.-6 ) .AND. ( IDIST.LT.1 .OR. IDIST.GT.4 ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLATM1', -INFO )
         RETURN
      END IF
*
      // Compute D according to COND and MODE
*
      IF( MODE.NE.0 ) THEN
         GO TO ( 10, 30, 50, 70, 90, 110 )ABS( MODE )
*
         // One large D value:
*
   10    CONTINUE
         DO 20 I = 1, N
            D( I ) = ONE / COND
   20    CONTINUE
         D( 1 ) = ONE
         GO TO 120
*
         // One small D value:
*
   30    CONTINUE
         DO 40 I = 1, N
            D( I ) = ONE
   40    CONTINUE
         D( N ) = ONE / COND
         GO TO 120
*
         // Exponentially distributed D values:
*
   50    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1 ) THEN
            ALPHA = COND**( -ONE / REAL( N-1 ) )
            DO 60 I = 2, N
               D( I ) = ALPHA**( I-1 )
   60       CONTINUE
         END IF
         GO TO 120
*
         // Arithmetically distributed D values:
*
   70    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1 ) THEN
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / REAL( N-1 )
            DO 80 I = 2, N
               D( I ) = REAL( N-I )*ALPHA + TEMP
   80       CONTINUE
         END IF
         GO TO 120
*
         // Randomly distributed D values on ( 1/COND , 1):
*
   90    CONTINUE
         ALPHA = LOG( ONE / COND )
         DO 100 I = 1, N
            D( I ) = EXP( ALPHA*SLARAN( ISEED ) )
  100    CONTINUE
         GO TO 120
*
         // Randomly distributed D values from IDIST
*
  110    CONTINUE
         CALL CLARNV( IDIST, ISEED, N, D )
*
  120    CONTINUE
*
         // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
         // random signs to D
*
         IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.1 ) THEN
            DO 130 I = 1, N
               CTEMP = CLARND( 3, ISEED )
               D( I ) = D( I )*( CTEMP / ABS( CTEMP ) )
  130       CONTINUE
         END IF
*
         // Reverse if MODE < 0
*
         IF( MODE.LT.0 ) THEN
            DO 140 I = 1, N / 2
               CTEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = CTEMP
  140       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
      // End of CLATM1
*
      END
