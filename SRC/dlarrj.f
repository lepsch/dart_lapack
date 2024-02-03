      SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST, RTOL, OFFSET, W, WERR, WORK, IWORK, PIVMIN, SPDIAM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IFIRST, ILAST, INFO, N, OFFSET;
      double             PIVMIN, RTOL, SPDIAM;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E2( * ), W( * ), WERR( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, HALF;
      const            ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0 ;
      int       MAXITR;
      // ..
      // .. Local Scalars ..
      int                CNT, I, I1, I2, II, ITER, J, K, NEXT, NINT, OLNINT, P, PREV, SAVI1;
      double             DPLUS, FAC, LEFT, MID, RIGHT, S, TMP, WIDTH;

      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if ( N.LE.0 ) {
         RETURN
      }

      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2

      // Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
      // The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
      // Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
      // for an unconverged interval is set to the index of the next unconverged
      // interval, and is -1 or 0 for a converged interval. Thus a linked
      // list of unconverged intervals is set up.


      I1 = IFIRST
      I2 = ILAST
      // The number of unconverged intervals
      NINT = 0
      // The last unconverged interval found
      PREV = 0
      for (I = I1; I <= I2; I++) { // 75
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         MID = W(II)
         RIGHT = W( II ) + WERR( II )
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )

         // The following test prevents the test of converged intervals
         if ( WIDTH.LT.RTOL*TMP ) {
            // This interval has already converged and does not need refinement.
            // (Note that the gaps might change through refining the
             // eigenvalues, however, they can only get bigger.)
            // Remove it from the list.
            IWORK( K-1 ) = -1
            // Make sure that I1 always points to the first unconverged interval
            IF((I.EQ.I1).AND.(I.LT.I2)) I1 = I + 1
            IF((PREV.GE.I1).AND.(I.LE.I2)) IWORK( 2*PREV-1 ) = I + 1
         } else {
            // unconverged interval found
            PREV = I
            // Make sure that [LEFT,RIGHT] contains the desired eigenvalue

            // Do while( CNT(LEFT).GT.I-1 )

            FAC = ONE
            } // 20
            CNT = 0
            S = LEFT
            DPLUS = D( 1 ) - S
            if (DPLUS.LT.ZERO) CNT = CNT + 1;
            for (J = 2; J <= N; J++) { // 30
               DPLUS = D( J ) - S - E2( J-1 )/DPLUS
               if (DPLUS.LT.ZERO) CNT = CNT + 1;
            } // 30
            if ( CNT.GT.I-1 ) {
               LEFT = LEFT - WERR( II )*FAC
               FAC = TWO*FAC
               GO TO 20
            }

            // Do while( CNT(RIGHT).LT.I )

            FAC = ONE
            } // 50
            CNT = 0
            S = RIGHT
            DPLUS = D( 1 ) - S
            if (DPLUS.LT.ZERO) CNT = CNT + 1;
            for (J = 2; J <= N; J++) { // 60
               DPLUS = D( J ) - S - E2( J-1 )/DPLUS
               if (DPLUS.LT.ZERO) CNT = CNT + 1;
            } // 60
            if ( CNT.LT.I ) {
               RIGHT = RIGHT + WERR( II )*FAC
               FAC = TWO*FAC
               GO TO 50
            }
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = CNT
         }
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
      } // 75


      SAVI1 = I1

      // Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
      // and while (ITER.LT.MAXITR)

      ITER = 0
      } // 80
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      for (P = 1; P <= OLNINT; P++) { // 100
         K = 2*I
         II = I - OFFSET
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT )

         // semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
          if ( ( WIDTH.LT.RTOL*TMP ) .OR. (ITER.EQ.MAXITR) ) {
            // reduce number of unconverged intervals
            NINT = NINT - 1
            // Mark interval as converged.
            IWORK( K-1 ) = 0
            if ( I1.EQ.I ) {
               I1 = NEXT
            } else {
               // Prev holds the last unconverged interval previously examined
               if (PREV.GE.I1) IWORK( 2*PREV-1 ) = NEXT;
            }
            I = NEXT
            GO TO 100
         }
         PREV = I

         // Perform one bisection step

         CNT = 0
         S = MID
         DPLUS = D( 1 ) - S
         if (DPLUS.LT.ZERO) CNT = CNT + 1;
         for (J = 2; J <= N; J++) { // 90
            DPLUS = D( J ) - S - E2( J-1 )/DPLUS
            if (DPLUS.LT.ZERO) CNT = CNT + 1;
         } // 90
         if ( CNT.LE.I-1 ) {
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
      IF( ( NINT.GT.0 ).AND.(ITER.LE.MAXITR) ) GO TO 80


      // At this point, all the intervals have converged
      for (I = SAVI1; I <= ILAST; I++) { // 110
         K = 2*I
         II = I - OFFSET
         // All intervals marked by '0' have been refined.
         if ( IWORK( K-1 ).EQ.0 ) {
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
         }
      } // 110


      RETURN

      // End of DLARRJ

      }
