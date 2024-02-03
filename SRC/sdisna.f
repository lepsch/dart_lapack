      SUBROUTINE SDISNA( JOB, M, N, D, SEP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOB;
      int                INFO, M, N;
*     ..
*     .. Array Arguments ..
      REAL               D( * ), SEP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               DECR, EIGEN, INCR, LEFT, RIGHT, SING;
      int                I, K;
      REAL               ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      EIGEN = LSAME( JOB, 'E' )
      LEFT = LSAME( JOB, 'L' )
      RIGHT = LSAME( JOB, 'R' )
      SING = LEFT .OR. RIGHT
      IF( EIGEN ) THEN
         K = M
      ELSE IF( SING ) THEN
         K = MIN( M, N )
      END IF
      IF( .NOT.EIGEN .AND. .NOT.SING ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( K.LT.0 ) THEN
         INFO = -3
      ELSE
         INCR = .TRUE.
         DECR = .TRUE.
         DO 10 I = 1, K - 1
            IF( INCR ) INCR = INCR .AND. D( I ).LE.D( I+1 )             IF( DECR ) DECR = DECR .AND. D( I ).GE.D( I+1 )
   10    CONTINUE
         IF( SING .AND. K.GT.0 ) THEN
            IF( INCR ) INCR = INCR .AND. ZERO.LE.D( 1 )             IF( DECR ) DECR = DECR .AND. D( K ).GE.ZERO
         END IF
         IF( .NOT.( INCR .OR. DECR ) ) INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SDISNA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 ) RETURN
*
*     Compute reciprocal condition numbers
*
      IF( K.EQ.1 ) THEN
         SEP( 1 ) = SLAMCH( 'O' )
      ELSE
         OLDGAP = ABS( D( 2 )-D( 1 ) )
         SEP( 1 ) = OLDGAP
         DO 20 I = 2, K - 1
            NEWGAP = ABS( D( I+1 )-D( I ) )
            SEP( I ) = MIN( OLDGAP, NEWGAP )
            OLDGAP = NEWGAP
   20    CONTINUE
         SEP( K ) = OLDGAP
      END IF
      IF( SING ) THEN
         IF( ( LEFT .AND. M.GT.N ) .OR. ( RIGHT .AND. M.LT.N ) ) THEN
            IF( INCR ) SEP( 1 ) = MIN( SEP( 1 ), D( 1 ) )             IF( DECR ) SEP( K ) = MIN( SEP( K ), D( K ) )
         END IF
      END IF
*
*     Ensure that reciprocal condition numbers are not less than
*     threshold, in order to limit the size of the error bound
*
      EPS = SLAMCH( 'E' )
      SAFMIN = SLAMCH( 'S' )
      ANORM = MAX( ABS( D( 1 ) ), ABS( D( K ) ) )
      IF( ANORM.EQ.ZERO ) THEN
         THRESH = EPS
      ELSE
         THRESH = MAX( EPS*ANORM, SAFMIN )
      END IF
      DO 30 I = 1, K
         SEP( I ) = MAX( SEP( I ), THRESH )
   30 CONTINUE
*
      RETURN
*
*     End of SDISNA
*
      END
