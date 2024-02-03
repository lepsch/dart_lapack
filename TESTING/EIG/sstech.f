      SUBROUTINE SSTECH( N, A, B, EIG, TOL, WORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      REAL               TOL
      // ..
      // .. Array Arguments ..
      REAL               A( * ), B( * ), EIG( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT;
      REAL               EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSTECT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Check input parameters

      INFO = 0
      if (N.EQ.0) RETURN;
      if ( N.LT.0 ) {
         INFO = -1
         RETURN
      }
      if ( TOL.LT.ZERO ) {
         INFO = -5
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      UNFLEP = SLAMCH( 'Safe minimum' ) / EPS
      EPS = TOL*EPS

      // Compute maximum absolute eigenvalue, error tolerance

      MX = ABS( EIG( 1 ) )
      for (I = 2; I <= N; I++) { // 10
         MX = MAX( MX, ABS( EIG( I ) ) )
      } // 10
      EPS = MAX( EPS*MX, UNFLEP )

      // Sort eigenvalues from EIG into WORK

      for (I = 1; I <= N; I++) { // 20
         WORK( I ) = EIG( I )
      } // 20
      for (I = 1; I <= N - 1; I++) { // 40
         ISUB = 1
         EMIN = WORK( 1 )
         for (J = 2; J <= N + 1 - I; J++) { // 30
            if ( WORK( J ).LT.EMIN ) {
               ISUB = J
               EMIN = WORK( J )
            }
         } // 30
         if ( ISUB.NE.N+1-I ) {
            WORK( ISUB ) = WORK( N+1-I )
            WORK( N+1-I ) = EMIN
         }
      } // 40

      // TPNT points to singular value at right endpoint of interval
      // BPNT points to singular value at left  endpoint of interval

      TPNT = 1
      BPNT = 1

      // Begin loop over all intervals

      } // 50
      UPPER = WORK( TPNT ) + EPS
      LOWER = WORK( BPNT ) - EPS

      // Begin loop merging overlapping intervals

      } // 60
      if (BPNT.EQ.N) GO TO 70;
      TUPPR = WORK( BPNT+1 ) + EPS
      if (TUPPR.LT.LOWER) GO TO 70;

      // Merge

      BPNT = BPNT + 1
      LOWER = WORK( BPNT ) - EPS
      GO TO 60
      } // 70

      // Count singular values in interval [ LOWER, UPPER ]

      sstect(N, A, B, LOWER, NUML );
      sstect(N, A, B, UPPER, NUMU );
      COUNT = NUMU - NUML
      if ( COUNT.NE.BPNT-TPNT+1 ) {

         // Wrong number of singular values in interval

         INFO = TPNT
         GO TO 80
      }
      TPNT = BPNT + 1
      BPNT = TPNT
      if (TPNT.LE.N) GO TO 50;
      } // 80
      RETURN

      // End of SSTECH

      }
