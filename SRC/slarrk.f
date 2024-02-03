      SUBROUTINE SLARRK( N, IW, GL, GU, D, E2, PIVMIN, RELTOL, W, WERR, INFO)
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int       INFO, IW, N;
      REAL                PIVMIN, RELTOL, GL, GU, W, WERR
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E2( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               FUDGE, HALF, TWO, ZERO
      PARAMETER          ( HALF = 0.5E0, TWO = 2.0E0, FUDGE = TWO, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      int       I, IT, ITMAX, NEGCNT;
      REAL               ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1, TMP2, TNORM
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL   SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         INFO = 0
         RETURN
      END IF
*
*     Get machine constants
      EPS = SLAMCH( 'P' )

      TNORM = MAX( ABS( GL ), ABS( GU ) )
      RTOLI = RELTOL
      ATOLI = FUDGE*TWO*PIVMIN
       ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2

      INFO = -1

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
      IT = 0

 10   CONTINUE
*
*     Check if interval converged or maximum number of iterations reached
*
      TMP1 = ABS( RIGHT - LEFT )
      TMP2 = MAX( ABS(RIGHT), ABS(LEFT) )
      IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) ) THEN
         INFO = 0
         GOTO 30
      ENDIF
      IF(IT.GT.ITMAX) GOTO 30

*
*     Count number of negative pivots for mid-point
*
      IT = IT + 1
      MID = HALF * (LEFT + RIGHT)
      NEGCNT = 0
      TMP1 = D( 1 ) - MID
      IF( ABS( TMP1 ).LT.PIVMIN ) TMP1 = -PIVMIN       IF( TMP1.LE.ZERO ) NEGCNT = NEGCNT + 1
*
      DO 20 I = 2, N
         TMP1 = D( I ) - E2( I-1 ) / TMP1 - MID
         IF( ABS( TMP1 ).LT.PIVMIN ) TMP1 = -PIVMIN          IF( TMP1.LE.ZERO ) NEGCNT = NEGCNT + 1
 20   CONTINUE

      IF(NEGCNT.GE.IW) THEN
         RIGHT = MID
      ELSE
         LEFT = MID
      ENDIF
      GOTO 10

 30   CONTINUE
*
*     Converged or maximum number of iterations reached
*
      W = HALF * (LEFT + RIGHT)
      WERR = HALF * ABS( RIGHT - LEFT )

      RETURN
*
*     End of SLARRK
*
      END
