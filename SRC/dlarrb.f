      SUBROUTINE DLARRB( N, D, LLD, IFIRST, ILAST, RTOL1, RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK, PIVMIN, SPDIAM, TWIST, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IFIRST, ILAST, INFO, N, OFFSET, TWIST;
      double             PIVMIN, RTOL1, RTOL2, SPDIAM;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), LLD( * ), W( * ), WERR( * ), WGAP( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, TWO, HALF;
      const            ZERO = 0.0D0, TWO = 2.0D0, HALF = 0.5D0 ;
      int       MAXITR;
      // ..
      // .. Local Scalars ..
      int                I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT, OLNINT, PREV, R;
      double             BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH, RGAP, RIGHT, TMP, WIDTH;
      // ..
      // .. External Functions ..
      int                DLANEG;
      // EXTERNAL DLANEG

      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if ( N <= 0 ) {
         RETURN
      }

      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
      MNWDTH = TWO * PIVMIN

      R = TWIST
      IF((R < 1) || (R > N)) R = N

      // Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
      // The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
      // Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
      // for an unconverged interval is set to the index of the next unconverged
      // interval, and is -1 or 0 for a converged interval. Thus a linked
      // list of unconverged intervals is set up.

      I1 = IFIRST
      // The number of unconverged intervals
      NINT = 0
      // The last unconverged interval found
      PREV = 0

      RGAP = WGAP( I1-OFFSET )
      for (I = I1; I <= ILAST; I++) { // 75
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         RIGHT = W( II ) + WERR( II )
         LGAP = RGAP
         RGAP = WGAP( II )
         GAP = MIN( LGAP, RGAP )

         // Make sure that [LEFT,RIGHT] contains the desired eigenvalue
         // Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT

         // Do while( NEGCNT(LEFT) > I-1 )

         BACK = WERR( II )
         } // 20
         NEGCNT = DLANEG( N, D, LLD, LEFT, PIVMIN, R )
         if ( NEGCNT > I-1 ) {
            LEFT = LEFT - BACK
            BACK = TWO*BACK
            GO TO 20
         }

         // Do while( NEGCNT(RIGHT) < I )
         // Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT

         BACK = WERR( II )
         } // 50

         NEGCNT = DLANEG( N, D, LLD, RIGHT, PIVMIN, R )
          if ( NEGCNT < I ) {
             RIGHT = RIGHT + BACK
             BACK = TWO*BACK
             GO TO 50
          }
         WIDTH = HALF*ABS( LEFT - RIGHT )
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         if ( WIDTH <= CVRGD || WIDTH <= MNWDTH ) {
            // This interval has already converged and does not need refinement.
            // (Note that the gaps might change through refining the
             // eigenvalues, however, they can only get bigger.)
            // Remove it from the list.
            IWORK( K-1 ) = -1
            // Make sure that I1 always points to the first unconverged interval
            IF((I == I1) && (I < ILAST)) I1 = I + 1
            IF((PREV >= I1) && (I <= ILAST)) IWORK( 2*PREV-1 ) = I + 1
         } else {
            // unconverged interval found
            PREV = I
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = NEGCNT
         }
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
      } // 75


      // Do while( NINT > 0 ), i.e. there are still unconverged intervals
      // and while (ITER < MAXITR)

      ITER = 0
      } // 80
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      for (IP = 1; IP <= OLNINT; IP++) { // 100
         K = 2*I
         II = I - OFFSET
         RGAP = WGAP( II )
         LGAP = RGAP
         if (II > 1) LGAP = WGAP( II-1 );
         GAP = MIN( LGAP, RGAP )
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT )

         // semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         if ( ( WIDTH <= CVRGD ) || ( WIDTH <= MNWDTH ) || ( ITER == MAXITR ) ) {
            // reduce number of unconverged intervals
            NINT = NINT - 1
            // Mark interval as converged.
            IWORK( K-1 ) = 0
            if ( I1 == I ) {
               I1 = NEXT
            } else {
               // Prev holds the last unconverged interval previously examined
               if (PREV >= I1) IWORK( 2*PREV-1 ) = NEXT;
            }
            I = NEXT
            GO TO 100
         }
         PREV = I

         // Perform one bisection step

         NEGCNT = DLANEG( N, D, LLD, MID, PIVMIN, R )
         if ( NEGCNT <= I-1 ) {
            WORK( K-1 ) = MID
         } else {
            WORK( K ) = MID
         }
         I = NEXT
      } // 100
      ITER = ITER + 1
      // do another loop if there are still unconverged intervals
      // However, in the last iteration, all intervals are accepted
      // since this is the best we can do.
      IF( ( NINT > 0 ) && (ITER <= MAXITR) ) GO TO 80


      // At this point, all the intervals have converged
      for (I = IFIRST; I <= ILAST; I++) { // 110
         K = 2*I
         II = I - OFFSET
         // All intervals marked by '0' have been refined.
         if ( IWORK( K-1 ) == 0 ) {
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
         }
      } // 110

      for (I = IFIRST+1; I <= ILAST; I++) { // 111
         K = 2*I
         II = I - OFFSET
         WGAP( II-1 ) = MAX( ZERO, W(II) - WERR (II) - W( II-1 ) - WERR( II-1 ))
      } // 111

      RETURN

      // End of DLARRB

      }
