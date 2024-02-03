      SUBROUTINE SLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, N, RANK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               COND
      int                IDIST, INFO, IRSIGN, MODE, N, RANK
*     ..
*     .. Array Arguments ..
      REAL               D( * )
      int                ISEED( 4 )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      REAL               ALPHA, TEMP
      int                I
*     ..
*     .. External Functions ..
      REAL               SLARAN
      EXTERNAL           SLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARNV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, EXP, LOG, REAL
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters. Initialize flags & seed.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Set INFO if an error
*
      IF( MODE.LT.-6 .OR. MODE.GT.6 ) THEN
         INFO = -1
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. ( IRSIGN.NE.0 .AND. IRSIGN.NE.1 ) ) THEN
         INFO = -2
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. COND.LT.ONE ) THEN
         INFO = -3
      ELSE IF( ( MODE.EQ.6 .OR. MODE.EQ.-6 ) .AND. ( IDIST.LT.1 .OR. IDIST.GT.3 ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLATM7', -INFO )
         RETURN
      END IF
*
*     Compute D according to COND and MODE
*
      IF( MODE.NE.0 ) THEN
         GO TO ( 100, 130, 160, 190, 210, 230 )ABS( MODE )
*
*        One large D value:
*
  100    CONTINUE
         DO 110 I = 2, RANK
            D( I ) = ONE / COND
  110    CONTINUE
         DO 120 I = RANK + 1, N
            D( I ) = ZERO
  120    CONTINUE
         D( 1 ) = ONE
         GO TO 240
*
*        One small D value:
*
  130    CONTINUE
         DO 140 I = 1, RANK - 1
            D( I ) = ONE
  140    CONTINUE
         DO 150 I = RANK + 1, N
            D( I ) = ZERO
  150    CONTINUE
         D( RANK ) = ONE / COND
         GO TO 240
*
*        Exponentially distributed D values:
*
  160    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1  .AND. RANK.GT.1 ) THEN
            ALPHA = COND**( -ONE / REAL( RANK-1 ) )
            DO 170 I = 2, RANK
               D( I ) = ALPHA**( I-1 )
  170       CONTINUE
            DO 180 I = RANK + 1, N
               D( I ) = ZERO
  180       CONTINUE
         END IF
         GO TO 240
*
*        Arithmetically distributed D values:
*
  190    CONTINUE
         D( 1 ) = ONE
         IF( N.GT.1 ) THEN
            TEMP = ONE / COND
            ALPHA = ( ONE-TEMP ) / REAL( N-1 )
            DO 200 I = 2, N
               D( I ) = REAL( N-I )*ALPHA + TEMP
  200       CONTINUE
         END IF
         GO TO 240
*
*        Randomly distributed D values on ( 1/COND , 1):
*
  210    CONTINUE
         ALPHA = LOG( ONE / COND )
         DO 220 I = 1, N
            D( I ) = EXP( ALPHA*SLARAN( ISEED ) )
  220    CONTINUE
         GO TO 240
*
*        Randomly distributed D values from IDIST
*
  230    CONTINUE
         CALL SLARNV( IDIST, ISEED, N, D )
*
  240    CONTINUE
*
*        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
*        random signs to D
*
         IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.1 ) THEN
            DO 250 I = 1, N
               TEMP = SLARAN( ISEED )
               IF( TEMP.GT.HALF ) D( I ) = -D( I )
  250       CONTINUE
         END IF
*
*        Reverse if MODE < 0
*
         IF( MODE.LT.0 ) THEN
            DO 260 I = 1, N / 2
               TEMP = D( I )
               D( I ) = D( N+1-I )
               D( N+1-I ) = TEMP
  260       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of SLATM7
*
      END
