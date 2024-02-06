      void sstech(N, A, B, EIG, TOL, WORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N;
      double               TOL;
      double               A( * ), B( * ), EIG( * ), WORK( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT;
      double               EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSTECT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      // Check input parameters

      INFO = 0;
      if (N == 0) return;
      if ( N < 0 ) {
         INFO = -1;
         return;
      }
      if ( TOL < ZERO ) {
         INFO = -5;
         return;
      }

      // Get machine constants

      EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      UNFLEP = SLAMCH( 'Safe minimum' ) / EPS;
      EPS = TOL*EPS;

      // Compute maximum absolute eigenvalue, error tolerance

      MX = ( EIG( 1 ) ).abs();
      for (I = 2; I <= N; I++) { // 10
         MX = max( MX, ( EIG( I ) ).abs() );
      } // 10
      EPS = max( EPS*MX, UNFLEP );

      // Sort eigenvalues from EIG into WORK

      for (I = 1; I <= N; I++) { // 20
         WORK[I] = EIG( I );
      } // 20
      for (I = 1; I <= N - 1; I++) { // 40
         ISUB = 1;
         EMIN = WORK( 1 );
         for (J = 2; J <= N + 1 - I; J++) { // 30
            if ( WORK( J ) < EMIN ) {
               ISUB = J;
               EMIN = WORK( J );
            }
         } // 30
         if ( ISUB != N+1-I ) {
            WORK[ISUB] = WORK( N+1-I );
            WORK[N+1-I] = EMIN;
         }
      } // 40

      // TPNT points to singular value at right endpoint of interval
      // BPNT points to singular value at left  endpoint of interval

      TPNT = 1;
      BPNT = 1;

      // Begin loop over all intervals

      } // 50
      UPPER = WORK( TPNT ) + EPS;
      LOWER = WORK( BPNT ) - EPS;

      // Begin loop merging overlapping intervals

      } // 60
      if (BPNT == N) GO TO 70;
      TUPPR = WORK( BPNT+1 ) + EPS;
      if (TUPPR < LOWER) GO TO 70;

      // Merge

      BPNT = BPNT + 1;
      LOWER = WORK( BPNT ) - EPS;
      GO TO 60;
      } // 70

      // Count singular values in interval [ LOWER, UPPER ]

      sstect(N, A, B, LOWER, NUML );
      sstect(N, A, B, UPPER, NUMU );
      COUNT = NUMU - NUML;
      if ( COUNT != BPNT-TPNT+1 ) {

         // Wrong number of singular values in interval

         INFO = TPNT;
         GO TO 80;
      }
      TPNT = BPNT + 1;
      BPNT = TPNT;
      if (TPNT <= N) GO TO 50;
      } // 80
      return;
      }
