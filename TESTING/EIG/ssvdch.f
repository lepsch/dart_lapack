      SUBROUTINE SSVDCH( N, S, E, SVD, TOL, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      REAL               TOL
      // ..
      // .. Array Arguments ..
      REAL               E( * ), S( * ), SVD( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                BPNT, COUNT, NUML, NUMU, TPNT;
      REAL               EPS, LOWER, OVFL, TUPPR, UNFL, UNFLEP, UPPER
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSVDCT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine constants

      INFO = 0
      IF( N.LE.0 ) RETURN
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      // UNFLEP is chosen so that when an eigenvalue is multiplied by the
      // scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in SSVDCT, it exceeds
      // sqrt(UNFL), which is the lower limit for SSVDCT.

      UNFLEP = ( SQRT( SQRT( UNFL ) ) / SQRT( OVFL ) )*SVD( 1 ) + UNFL / EPS

      // The value of EPS works best when TOL .GE. 10.

      EPS = TOL*MAX( N / 10, 1 )*EPS

      // TPNT points to singular value at right endpoint of interval
      // BPNT points to singular value at left  endpoint of interval

      TPNT = 1
      BPNT = 1

      // Begin loop over all intervals

      } // 10
      UPPER = ( ONE+EPS )*SVD( TPNT ) + UNFLEP
      LOWER = ( ONE-EPS )*SVD( BPNT ) - UNFLEP
      IF( LOWER.LE.UNFLEP ) LOWER = -UPPER

      // Begin loop merging overlapping intervals

      } // 20
      IF( BPNT.EQ.N ) GO TO 30
      TUPPR = ( ONE+EPS )*SVD( BPNT+1 ) + UNFLEP
      IF( TUPPR.LT.LOWER ) GO TO 30

      // Merge

      BPNT = BPNT + 1
      LOWER = ( ONE-EPS )*SVD( BPNT ) - UNFLEP
      IF( LOWER.LE.UNFLEP ) LOWER = -UPPER
      GO TO 20
      } // 30

      // Count singular values in interval [ LOWER, UPPER ]

      ssvdct(N, S, E, LOWER, NUML );
      ssvdct(N, S, E, UPPER, NUMU );
      COUNT = NUMU - NUML
      IF( LOWER.LT.ZERO ) COUNT = COUNT / 2
      if ( COUNT.NE.BPNT-TPNT+1 ) {

         // Wrong number of singular values in interval

         INFO = TPNT
         GO TO 40
      }
      TPNT = BPNT + 1
      BPNT = TPNT
      IF( TPNT.LE.N ) GO TO 10
      } // 40
      RETURN

      // End of SSVDCH

      }
