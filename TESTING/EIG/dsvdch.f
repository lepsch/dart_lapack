      SUBROUTINE DSVDCH( N, S, E, SVD, TOL, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             TOL;
      // ..
      // .. Array Arguments ..
      double             E( * ), S( * ), SVD( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D0 ;
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                BPNT, COUNT, NUML, NUMU, TPNT;
      double             EPS, LOWER, OVFL, TUPPR, UNFL, UNFLEP, UPPER;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSVDCT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine constants

      INFO = 0
      IF( N.LE.0 ) RETURN
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      EPS = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // UNFLEP is chosen so that when an eigenvalue is multiplied by the
      // scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in DSVDCT, it exceeds
      // sqrt(UNFL), which is the lower limit for DSVDCT.

      UNFLEP = ( SQRT( SQRT( UNFL ) ) / SQRT( OVFL ) )*SVD( 1 ) + UNFL / EPS

      // The value of EPS works best when TOL .GE. 10.

      EPS = TOL*MAX( N / 10, 1 )*EPS

      // TPNT points to singular value at right endpoint of interval
      // BPNT points to singular value at left  endpoint of interval

      TPNT = 1
      BPNT = 1

      // Begin loop over all intervals

   10 CONTINUE
      UPPER = ( ONE+EPS )*SVD( TPNT ) + UNFLEP
      LOWER = ( ONE-EPS )*SVD( BPNT ) - UNFLEP
      IF( LOWER.LE.UNFLEP ) LOWER = -UPPER

      // Begin loop merging overlapping intervals

   20 CONTINUE
      IF( BPNT.EQ.N ) GO TO 30
      TUPPR = ( ONE+EPS )*SVD( BPNT+1 ) + UNFLEP
      IF( TUPPR.LT.LOWER ) GO TO 30

      // Merge

      BPNT = BPNT + 1
      LOWER = ( ONE-EPS )*SVD( BPNT ) - UNFLEP
      IF( LOWER.LE.UNFLEP ) LOWER = -UPPER
      GO TO 20
   30 CONTINUE

      // Count singular values in interval [ LOWER, UPPER ]

      CALL DSVDCT( N, S, E, LOWER, NUML )
      CALL DSVDCT( N, S, E, UPPER, NUMU )
      COUNT = NUMU - NUML
      IF( LOWER.LT.ZERO ) COUNT = COUNT / 2
      IF( COUNT.NE.BPNT-TPNT+1 ) THEN

         // Wrong number of singular values in interval

         INFO = TPNT
         GO TO 40
      END IF
      TPNT = BPNT + 1
      BPNT = TPNT
      IF( TPNT.LE.N ) GO TO 10
   40 CONTINUE
      RETURN

      // End of DSVDCH

      }
