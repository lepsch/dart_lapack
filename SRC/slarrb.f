      SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1, RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK, PIVMIN, SPDIAM, TWIST, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IFIRST, ILAST, INFO, N, OFFSET, TWIST
      REAL               PIVMIN, RTOL1, RTOL2, SPDIAM
*     ..
*     .. Array Arguments ..
      int                IWORK( * )
      REAL               D( * ), LLD( * ), W( * ), WERR( * ), WGAP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, TWO, HALF
      PARAMETER        ( ZERO = 0.0E0, TWO = 2.0E0, HALF = 0.5E0 )
      int       MAXITR
*     ..
*     .. Local Scalars ..
      int                I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT, OLNINT, PREV, R       REAL               BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH, RGAP, RIGHT, TMP, WIDTH
*     ..
*     .. External Functions ..
      int                SLANEG
      EXTERNAL           SLANEG
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RETURN
      END IF
*
      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
      MNWDTH = TWO * PIVMIN
*
      R = TWIST
      IF((R.LT.1).OR.(R.GT.N)) R = N
*
*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
*     for an unconverged interval is set to the index of the next unconverged
*     interval, and is -1 or 0 for a converged interval. Thus a linked
*     list of unconverged intervals is set up.
*
      I1 = IFIRST
*     The number of unconverged intervals
      NINT = 0
*     The last unconverged interval found
      PREV = 0

      RGAP = WGAP( I1-OFFSET )
      DO 75 I = I1, ILAST
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         RIGHT = W( II ) + WERR( II )
         LGAP = RGAP
         RGAP = WGAP( II )
         GAP = MIN( LGAP, RGAP )

*        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
*
*        Do while( NEGCNT(LEFT).GT.I-1 )
*
         BACK = WERR( II )
 20      CONTINUE
         NEGCNT = SLANEG( N, D, LLD, LEFT, PIVMIN, R )
         IF( NEGCNT.GT.I-1 ) THEN
            LEFT = LEFT - BACK
            BACK = TWO*BACK
            GO TO 20
         END IF
*
*        Do while( NEGCNT(RIGHT).LT.I )
*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
*
         BACK = WERR( II )
 50      CONTINUE

         NEGCNT = SLANEG( N, D, LLD, RIGHT, PIVMIN, R )
          IF( NEGCNT.LT.I ) THEN
             RIGHT = RIGHT + BACK
             BACK = TWO*BACK
             GO TO 50
          END IF
         WIDTH = HALF*ABS( LEFT - RIGHT )
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( WIDTH.LE.CVRGD .OR. WIDTH.LE.MNWDTH ) THEN
*           This interval has already converged and does not need refinement.
*           (Note that the gaps might change through refining the
*            eigenvalues, however, they can only get bigger.)
*           Remove it from the list.
            IWORK( K-1 ) = -1
*           Make sure that I1 always points to the first unconverged interval
            IF((I.EQ.I1).AND.(I.LT.ILAST)) I1 = I + 1
            IF((PREV.GE.I1).AND.(I.LE.ILAST)) IWORK( 2*PREV-1 ) = I + 1
         ELSE
*           unconverged interval found
            PREV = I
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = NEGCNT
         END IF
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
 75   CONTINUE

*
*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
*     and while (ITER.LT.MAXITR)
*
      ITER = 0
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 IP = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         RGAP = WGAP( II )
         LGAP = RGAP
         IF(II.GT.1) LGAP = WGAP( II-1 )
         GAP = MIN( LGAP, RGAP )
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT )

*        semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( ( WIDTH.LE.CVRGD ) .OR. ( WIDTH.LE.MNWDTH ).OR. ( ITER.EQ.MAXITR ) )THEN
*           reduce number of unconverged intervals
            NINT = NINT - 1
*           Mark interval as converged.
            IWORK( K-1 ) = 0
            IF( I1.EQ.I ) THEN
               I1 = NEXT
            ELSE
*              Prev holds the last unconverged interval previously examined
               IF(PREV.GE.I1) IWORK( 2*PREV-1 ) = NEXT
            END IF
            I = NEXT
            GO TO 100
         END IF
         PREV = I
*
*        Perform one bisection step
*
         NEGCNT = SLANEG( N, D, LLD, MID, PIVMIN, R )
         IF( NEGCNT.LE.I-1 ) THEN
            WORK( K-1 ) = MID
         ELSE
            WORK( K ) = MID
         END IF
         I = NEXT
 100  CONTINUE
      ITER = ITER + 1
*     do another loop if there are still unconverged intervals
*     However, in the last iteration, all intervals are accepted
*     since this is the best we can do.
      IF( ( NINT.GT.0 ).AND.(ITER.LE.MAXITR) ) GO TO 80
*
*
*     At this point, all the intervals have converged
      DO 110 I = IFIRST, ILAST
         K = 2*I
         II = I - OFFSET
*        All intervals marked by '0' have been refined.
         IF( IWORK( K-1 ).EQ.0 ) THEN
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
         END IF
 110  CONTINUE
*
      DO 111 I = IFIRST+1, ILAST
         K = 2*I
         II = I - OFFSET
         WGAP( II-1 ) = MAX( ZERO, W(II) - WERR (II) - W( II-1 ) - WERR( II-1 ))
 111  CONTINUE

      RETURN
*
*     End of SLARRB
*
      END
