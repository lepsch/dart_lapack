      void dbdsqr(UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      double             NEGONE;
      const              NEGONE = -1.0 ;
      double             HNDRTH;
      const              HNDRTH = 0.01 ;
      double             TEN;
      const              TEN = 10.0 ;
      double             HNDRD;
      const              HNDRD = 100.0 ;
      double             MEIGTH;
      const              MEIGTH = -0.125 ;
      int                MAXITR;
      const              MAXITR = 6 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, ROTATE;
      int                I, IDIR, ISUB, ITER, ITERDIVN, J, LL, LLL, M, MAXITDIVN, NM1, NM12, NM13, OLDLL, OLDM;
      double             ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, SINR, SLL, SMAX, SMIN, SMINOA, SN, THRESH, TOL, TOLMUL, UNFL;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT, DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      LOWER = LSAME( UPLO, 'L' );
      if ( !LSAME( UPLO, 'U' ) && !LOWER ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NCVT < 0 ) {
         INFO = -3;
      } else if ( NRU < 0 ) {
         INFO = -4;
      } else if ( NCC < 0 ) {
         INFO = -5;
      } else if ( ( NCVT == 0 && LDVT < 1 ) || ( NCVT > 0 && LDVT < max( 1, N ) ) ) {
         INFO = -9;
      } else if ( LDU < max( 1, NRU ) ) {
         INFO = -11;
      } else if ( ( NCC == 0 && LDC < 1 ) || ( NCC > 0 && LDC < max( 1, N ) ) ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('DBDSQR', -INFO );
         return;
      }
      if (N == 0) return;
      IF( N == 1 ) GO TO 160;

      // ROTATE is true if any singular vectors desired, false otherwise

      ROTATE = ( NCVT > 0 ) || ( NRU > 0 ) || ( NCC > 0 );

      // If no singular vectors desired, use qd algorithm

      if ( !ROTATE ) {
         dlasq1(N, D, E, WORK, INFO );

      // If INFO equals 2, dqds didn't finish, try to finish

         if (INFO != 2) return;
         INFO = 0;
      }

      NM1 = N - 1;
      NM12 = NM1 + NM1;
      NM13 = NM12 + NM1;
      IDIR = 0;

      // Get machine constants

      EPS = DLAMCH( 'Epsilon' );
      UNFL = DLAMCH( 'Safe minimum' );

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      if ( LOWER ) {
         for (I = 1; I <= N - 1; I++) { // 10
            dlartg(D( I ), E( I ), CS, SN, R );
            D[I] = R;
            E[I] = SN*D( I+1 );
            D[I+1] = CS*D( I+1 );
            WORK[I] = CS;
            WORK[NM1+I] = SN;
         } // 10

         // Update singular vectors if desired

         if (NRU > 0) dlasr( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U, LDU );
         IF( NCC > 0 ) dlasr( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C, LDC );
      }

      // Compute singular values to relative accuracy TOL
      // (By setting TOL to be negative, algorithm will compute
      // singular values to absolute accuracy ABS(TOL)*norm(input matrix))

      TOLMUL = max( TEN, min( HNDRD, EPS**MEIGTH ) );
      TOL = TOLMUL*EPS;

      // Compute approximate maximum, minimum singular values

      SMAX = ZERO;
      for (I = 1; I <= N; I++) { // 20
         SMAX = max( SMAX, ( D( I ) ) ).abs();
      } // 20
      for (I = 1; I <= N - 1; I++) { // 30
         SMAX = max( SMAX, ( E( I ) ) ).abs();
      } // 30
      SMIN = ZERO;
      if ( TOL >= ZERO ) {

         // Relative accuracy desired

         SMINOA = ( D( 1 ) ).abs();
         if (SMINOA == ZERO) GO TO 50;
         MU = SMINOA;
         for (I = 2; I <= N; I++) { // 40
            MU = ( D( I ) ).abs()*( MU / ( MU+( E( I-1 ) ) ) ).abs();
            SMINOA = min( SMINOA, MU );
            if (SMINOA == ZERO) GO TO 50;
         } // 40
         } // 50
         SMINOA = SMINOA / sqrt( N.toDouble() );
         THRESH = max( TOL*SMINOA, MAXITR*(N*(N*UNFL)) );
      } else {

         // Absolute accuracy desired

         THRESH = max( ( TOL ).abs()*SMAX, MAXITR*(N*(N*UNFL)) );
      }

      // Prepare for main iteration loop for the singular values
      // (MAXIT is the maximum number of passes through the inner
      // loop permitted before nonconvergence signalled.)

      MAXITDIVN = MAXITR*N;
      ITERDIVN = 0;
      ITER = -1;
      OLDLL = -1;
      OLDM = -1;

      // M points to last element of unconverged part of matrix

      M = N;

      // Begin main iteration loop

      } // 60

      // Check for convergence or exceeding iteration count

      if (M <= 1) GO TO 160;

      if ( ITER >= N ) {
         ITER = ITER - N;
         ITERDIVN = ITERDIVN + 1;
         if (ITERDIVN >= MAXITDIVN) GO TO 200;
      }

      // Find diagonal block of matrix to work on

      if( TOL < ZERO && ( D( M ) ).abs() <= THRESH ) D( M ) = ZERO;
      SMAX = ( D( M ) ).abs();
      for (LLL = 1; LLL <= M - 1; LLL++) { // 70
         LL = M - LLL;
         ABSS = ( D( LL ) ).abs();
         ABSE = ( E( LL ) ).abs();
         if (TOL < ZERO && ABSS <= THRESH) D( LL ) = ZERO;
         IF( ABSE <= THRESH ) GO TO 80;
         SMAX = max( SMAX, ABSS, ABSE );
      } // 70
      LL = 0;
      GO TO 90;
      } // 80
      E[LL] = ZERO;

      // Matrix splits since E(LL) = 0

      if ( LL == M-1 ) {

         // Convergence of bottom singular value, return to top of loop

         M = M - 1;
         GO TO 60;
      }
      } // 90
      LL = LL + 1;

      // E(LL) through E(M-1) are nonzero, E(LL-1) is zero

      if ( LL == M-1 ) {

         // 2 by 2 block, handle separately

         dlasv2(D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, COSR, SINL, COSL );
         D[M-1] = SIGMX;
         E[M-1] = ZERO;
         D[M] = SIGMN;

         // Compute singular vectors, if desired

         if (NCVT > 0) drot( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, SINR );
         if( NRU > 0 ) drot( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL );
         IF( NCC > 0 ) drot( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, SINL );
         M = M - 2;
         GO TO 60;
      }

      // If working on new submatrix, choose shift direction
      // (from larger end diagonal element towards smaller)

      if ( LL > OLDM || M < OLDLL ) {
         if ( ( D( LL ) ).abs() >= ( D( M ) ) ).abs() {

            // Chase bulge from top (big end) to bottom (small end)

            IDIR = 1;
         } else {

            // Chase bulge from bottom (big end) to top (small end)

            IDIR = 2;
         }
      }

      // Apply convergence tests

      if ( IDIR == 1 ) {

         // Run convergence test in forward direction
         // First apply standard test to bottom of matrix

         if ( ( E( M-1 ) ).abs() <= ( TOL ).abs()*( D( M ) ).abs() || ( TOL < ZERO && ( E( M-1 ) ) <= THRESH ) ).abs() {
            E[M-1] = ZERO;
            GO TO 60;
         }

         if ( TOL >= ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion forward

            MU = ( D( LL ) ).abs();
            SMIN = MU;
            for (LLL = LL; LLL <= M - 1; LLL++) { // 100
               if ( ( E( LLL ) ).abs() <= TOL*MU ) {
                  E[LLL] = ZERO;
                  GO TO 60;
               }
               MU = ( D( LLL+1 ) ).abs()*( MU / ( MU+( E( LLL ) ) ) ).abs();
               SMIN = min( SMIN, MU );
            } // 100
         }

      } else {

         // Run convergence test in backward direction
         // First apply standard test to top of matrix

         if ( ( E( LL ) ).abs() <= ( TOL ).abs()*( D( LL ) ).abs() || ( TOL < ZERO && ( E( LL ) ) <= THRESH ) ).abs() {
            E[LL] = ZERO;
            GO TO 60;
         }

         if ( TOL >= ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion backward

            MU = ( D( M ) ).abs();
            SMIN = MU;
            for (LLL = M - 1; LLL >= LL; LLL--) { // 110
               if ( ( E( LLL ) ).abs() <= TOL*MU ) {
                  E[LLL] = ZERO;
                  GO TO 60;
               }
               MU = ( D( LLL ) ).abs()*( MU / ( MU+( E( LLL ) ) ) ).abs();
               SMIN = min( SMIN, MU );
            } // 110
         }
      }
      OLDLL = LL;
      OLDM = M;

      // Compute shift.  First, test if shifting would ruin relative
      // accuracy, and if so set the shift to zero.

      if ( TOL >= ZERO && N*TOL*( SMIN / SMAX ) <= max( EPS, HNDRTH*TOL ) ) {

         // Use a zero shift to avoid loss of relative accuracy

         SHIFT = ZERO;
      } else {

         // Compute the shift from 2-by-2 block at end of matrix

         if ( IDIR == 1 ) {
            SLL = ( D( LL ) ).abs();
            dlas2(D( M-1 ), E( M-1 ), D( M ), SHIFT, R );
         } else {
            SLL = ( D( M ) ).abs();
            dlas2(D( LL ), E( LL ), D( LL+1 ), SHIFT, R );
         }

         // Test if shift negligible, and if so set to zero

         if ( SLL > ZERO ) {
            if( ( SHIFT / SLL )**2 < EPS ) SHIFT = ZERO;
         }
      }

      // Increment iteration count

      ITER = ITER + M - LL;

      // If SHIFT = 0, do simplified QR iteration

      if ( SHIFT == ZERO ) {
         if ( IDIR == 1 ) {

            // Chase bulge from top to bottom
            // Save cosines and sines for later singular vector updates

            CS = ONE;
            OLDCS = ONE;
            for (I = LL; I <= M - 1; I++) { // 120
               dlartg(D( I )*CS, E( I ), CS, SN, R );
               if (I > LL) E( I-1 ) = OLDSN*R;
               dlartg(OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) );
               WORK[I-LL+1] = CS;
               WORK[I-LL+1+NM1] = SN;
               WORK[I-LL+1+NM12] = OLDCS;
               WORK[I-LL+1+NM13] = OLDSN;
            } // 120
            H = D( M )*CS;
            D[M] = H*OLDCS;
            E[M-1] = H*OLDSN;

            // Update singular vectors

            if (NCVT > 0) dlasr( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), WORK( N ), VT( LL, 1 ), LDVT );
            if( NRU > 0 ) dlasr( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), WORK( NM13+1 ), U( 1, LL ), LDU );
            IF( NCC > 0 ) dlasr( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), WORK( NM13+1 ), C( LL, 1 ), LDC );

            // Test convergence

            if( ( E( M-1 ) ).abs() <= THRESH ) E( M-1 ) = ZERO;

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            CS = ONE;
            OLDCS = ONE;
            for (I = M; I >= LL + 1; I--) { // 130
               dlartg(D( I )*CS, E( I-1 ), CS, SN, R );
               if (I < M) E( I ) = OLDSN*R;
               dlartg(OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) );
               WORK[I-LL] = CS;
               WORK[I-LL+NM1] = -SN;
               WORK[I-LL+NM12] = OLDCS;
               WORK[I-LL+NM13] = -OLDSN;
            } // 130
            H = D( LL )*CS;
            D[LL] = H*OLDCS;
            E[LL] = H*OLDSN;

            // Update singular vectors

            if (NCVT > 0) dlasr( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), WORK( NM13+1 ), VT( LL, 1 ), LDVT );
            if( NRU > 0 ) dlasr( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), WORK( N ), U( 1, LL ), LDU );
            IF( NCC > 0 ) dlasr( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), WORK( N ), C( LL, 1 ), LDC );

            // Test convergence

            if( ( E( LL ) ).abs() <= THRESH ) E( LL ) = ZERO;
         }
      } else {

         // Use nonzero shift

         if ( IDIR == 1 ) {

            // Chase bulge from top to bottom
            // Save cosines and sines for later singular vector updates

            F = ( ( D( LL ) ).abs()-SHIFT )* ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) );
            G = E( LL );
            for (I = LL; I <= M - 1; I++) { // 140
               dlartg(F, G, COSR, SINR, R );
               if (I > LL) E( I-1 ) = R;
               F = COSR*D( I ) + SINR*E( I );
               E[I] = COSR*E( I ) - SINR*D( I );
               G = SINR*D( I+1 );
               D[I+1] = COSR*D( I+1 );
               dlartg(F, G, COSL, SINL, R );
               D[I] = R;
               F = COSL*E( I ) + SINL*D( I+1 );
               D[I+1] = COSL*D( I+1 ) - SINL*E( I );
               if ( I < M-1 ) {
                  G = SINL*E( I+1 );
                  E[I+1] = COSL*E( I+1 );
               }
               WORK[I-LL+1] = COSR;
               WORK[I-LL+1+NM1] = SINR;
               WORK[I-LL+1+NM12] = COSL;
               WORK[I-LL+1+NM13] = SINL;
            } // 140
            E[M-1] = F;

            // Update singular vectors

            if (NCVT > 0) dlasr( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), WORK( N ), VT( LL, 1 ), LDVT );
            if( NRU > 0 ) dlasr( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), WORK( NM13+1 ), U( 1, LL ), LDU );
            IF( NCC > 0 ) dlasr( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), WORK( NM13+1 ), C( LL, 1 ), LDC );

            // Test convergence

            if( ( E( M-1 ) ).abs() <= THRESH ) E( M-1 ) = ZERO;

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            F = ( ( D( M ) ).abs()-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / D( M ) );
            G = E( M-1 );
            for (I = M; I >= LL + 1; I--) { // 150
               dlartg(F, G, COSR, SINR, R );
               if (I < M) E( I ) = R;
               F = COSR*D( I ) + SINR*E( I-1 );
               E[I-1] = COSR*E( I-1 ) - SINR*D( I );
               G = SINR*D( I-1 );
               D[I-1] = COSR*D( I-1 );
               dlartg(F, G, COSL, SINL, R );
               D[I] = R;
               F = COSL*E( I-1 ) + SINL*D( I-1 );
               D[I-1] = COSL*D( I-1 ) - SINL*E( I-1 );
               if ( I > LL+1 ) {
                  G = SINL*E( I-2 );
                  E[I-2] = COSL*E( I-2 );
               }
               WORK[I-LL] = COSR;
               WORK[I-LL+NM1] = -SINR;
               WORK[I-LL+NM12] = COSL;
               WORK[I-LL+NM13] = -SINL;
            } // 150
            E[LL] = F;

            // Test convergence

            if( ( E( LL ) ).abs() <= THRESH ) E( LL ) = ZERO;

            // Update singular vectors if desired

            if (NCVT > 0) dlasr( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), WORK( NM13+1 ), VT( LL, 1 ), LDVT );
            if( NRU > 0 ) dlasr( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), WORK( N ), U( 1, LL ), LDU );
            IF( NCC > 0 ) dlasr( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), WORK( N ), C( LL, 1 ), LDC );
         }
      }

      // QR iteration finished, go back and check convergence

      GO TO 60;

      // All singular values converged, so make them positive

      } // 160
      for (I = 1; I <= N; I++) { // 170
         if ( D( I ) < ZERO ) {
            D[I] = -D( I );

            // Change sign of singular vectors, if desired

            if (NCVT > 0) dscal( NCVT, NEGONE, VT( I, 1 ), LDVT );
         }
      } // 170

      // Sort the singular values into decreasing order (insertion sort on
      // singular values, but only one transposition per singular vector)

      for (I = 1; I <= N - 1; I++) { // 190

         // Scan for smallest D(I)

         ISUB = 1;
         SMIN = D( 1 );
         for (J = 2; J <= N + 1 - I; J++) { // 180
            if ( D( J ) <= SMIN ) {
               ISUB = J;
               SMIN = D( J );
            }
         } // 180
         if ( ISUB != N+1-I ) {

            // Swap singular values and vectors

            D[ISUB] = D( N+1-I );
            D[N+1-I] = SMIN;
            if (NCVT > 0) dswap( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), LDVT );
            if (NRU > 0) dswap( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 );
            IF( NCC > 0 ) dswap( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC );
         }
      } // 190
      GO TO 220;

      // Maximum number of iterations exceeded, failure to converge

      } // 200
      INFO = 0;
      for (I = 1; I <= N - 1; I++) { // 210
         if( E( I ) != ZERO ) INFO = INFO + 1;
      } // 210
      } // 220
      return;
      }
