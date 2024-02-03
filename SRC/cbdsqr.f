      SUBROUTINE CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, RWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RWORK( * );
      COMPLEX            C( LDC, * ), U( LDU, * ), VT( LDVT, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      REAL               ONE;
      const              ONE = 1.0 ;
      REAL               NEGONE;
      const              NEGONE = -1.0 ;
      REAL               HNDRTH;
      const              HNDRTH = 0.01 ;
      REAL               TEN;
      const              TEN = 10.0 ;
      REAL               HNDRD;
      const              HNDRD = 100.0 ;
      REAL               MEIGTH;
      const              MEIGTH = -0.125 ;
      int                MAXITR;
      const              MAXITR = 6 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, ROTATE;
      int                I, IDIR, ISUB, ITER, ITERDIVN, J, LL, LLL, M, MAXITDIVN, NM1, NM12, NM13, OLDLL, OLDM;
      REAL               ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, SINR, SLL, SMAX, SMIN, SMINOA, SN, THRESH, TOL, TOLMUL, UNFL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASR, CSROT, CSSCAL, CSWAP, SLARTG, SLAS2, SLASQ1, SLASV2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN, SQRT
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
      } else if ( ( NCVT == 0 && LDVT < 1 ) || ( NCVT > 0 && LDVT < MAX( 1, N ) ) ) {
         INFO = -9;
      } else if ( LDU < MAX( 1, NRU ) ) {
         INFO = -11;
      } else if ( ( NCC == 0 && LDC < 1 ) || ( NCC > 0 && LDC < MAX( 1, N ) ) ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('CBDSQR', -INFO );
         return;
      }
      if (N == 0) RETURN       IF( N == 1 ) GO TO 160;

      // ROTATE is true if any singular vectors desired, false otherwise

      ROTATE = ( NCVT > 0 ) || ( NRU > 0 ) || ( NCC > 0 );

      // If no singular vectors desired, use qd algorithm

      if ( !ROTATE ) {
         slasq1(N, D, E, RWORK, INFO );

      // If INFO equals 2, dqds didn't finish, try to finish

         if (INFO != 2) RETURN;
         INFO = 0;
      }

      NM1 = N - 1;
      NM12 = NM1 + NM1;
      NM13 = NM12 + NM1;
      IDIR = 0;

      // Get machine constants

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      if ( LOWER ) {
         for (I = 1; I <= N - 1; I++) { // 10
            slartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R;
            E( I ) = SN*D( I+1 );
            D( I+1 ) = CS*D( I+1 );
            RWORK( I ) = CS;
            RWORK( NM1+I ) = SN;
         } // 10

         // Update singular vectors if desired

         if (NRU > 0) CALL CLASR( 'R', 'V', 'F', NRU, N, RWORK( 1 ), RWORK( N ), U, LDU )          IF( NCC > 0 ) CALL CLASR( 'L', 'V', 'F', N, NCC, RWORK( 1 ), RWORK( N ), C, LDC );
      }

      // Compute singular values to relative accuracy TOL
      // (By setting TOL to be negative, algorithm will compute
      // singular values to absolute accuracy ABS(TOL)*norm(input matrix))

      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) );
      TOL = TOLMUL*EPS;

      // Compute approximate maximum, minimum singular values

      SMAX = ZERO;
      for (I = 1; I <= N; I++) { // 20
         SMAX = MAX( SMAX, ABS( D( I ) ) );
      } // 20
      for (I = 1; I <= N - 1; I++) { // 30
         SMAX = MAX( SMAX, ABS( E( I ) ) );
      } // 30
      SMIN = ZERO;
      if ( TOL >= ZERO ) {

         // Relative accuracy desired

         SMINOA = ABS( D( 1 ) );
         if (SMINOA == ZERO) GO TO 50;
         MU = SMINOA;
         for (I = 2; I <= N; I++) { // 40
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) );
            SMINOA = MIN( SMINOA, MU );
            if (SMINOA == ZERO) GO TO 50;
         } // 40
         } // 50
         SMINOA = SMINOA / SQRT( REAL( N ) );
         THRESH = MAX( TOL*SMINOA, MAXITR*(N*(N*UNFL)) );
      } else {

         // Absolute accuracy desired

         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*(N*(N*UNFL)) );
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

      IF( TOL < ZERO && ABS( D( M ) ) <= THRESH ) D( M ) = ZERO;
      SMAX = ABS( D( M ) );
      for (LLL = 1; LLL <= M - 1; LLL++) { // 70
         LL = M - LLL;
         ABSS = ABS( D( LL ) );
         ABSE = ABS( E( LL ) );
         if (TOL < ZERO && ABSS <= THRESH) D( LL ) = ZERO          IF( ABSE <= THRESH ) GO TO 80;
         SMAX = MAX( SMAX, ABSS, ABSE );
      } // 70
      LL = 0;
      GO TO 90;
      } // 80
      E( LL ) = ZERO;

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

         slasv2(D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, COSR, SINL, COSL );
         D( M-1 ) = SIGMX;
         E( M-1 ) = ZERO;
         D( M ) = SIGMN;

         // Compute singular vectors, if desired

         if (NCVT > 0) CALL CSROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, SINR )          IF( NRU > 0 ) CALL CSROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )          IF( NCC > 0 ) CALL CSROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, SINL );
         M = M - 2;
         GO TO 60;
      }

      // If working on new submatrix, choose shift direction
      // (from larger end diagonal element towards smaller)

      if ( LL > OLDM || M < OLDLL ) {
         if ( ABS( D( LL ) ) >= ABS( D( M ) ) ) {

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

         if ( ABS( E( M-1 ) ) <= ABS( TOL )*ABS( D( M ) ) || ( TOL < ZERO && ABS( E( M-1 ) ) <= THRESH ) ) {
            E( M-1 ) = ZERO;
            GO TO 60;
         }

         if ( TOL >= ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion forward

            MU = ABS( D( LL ) );
            SMIN = MU;
            for (LLL = LL; LLL <= M - 1; LLL++) { // 100
               if ( ABS( E( LLL ) ) <= TOL*MU ) {
                  E( LLL ) = ZERO;
                  GO TO 60;
               }
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) );
               SMIN = MIN( SMIN, MU );
            } // 100
         }

      } else {

         // Run convergence test in backward direction
         // First apply standard test to top of matrix

         if ( ABS( E( LL ) ) <= ABS( TOL )*ABS( D( LL ) ) || ( TOL < ZERO && ABS( E( LL ) ) <= THRESH ) ) {
            E( LL ) = ZERO;
            GO TO 60;
         }

         if ( TOL >= ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion backward

            MU = ABS( D( M ) );
            SMIN = MU;
            DO 110 LLL = M - 1, LL, -1;
               if ( ABS( E( LLL ) ) <= TOL*MU ) {
                  E( LLL ) = ZERO;
                  GO TO 60;
               }
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) );
               SMIN = MIN( SMIN, MU );
            } // 110
         }
      }
      OLDLL = LL;
      OLDM = M;

      // Compute shift.  First, test if shifting would ruin relative
      // accuracy, and if so set the shift to zero.

      if ( TOL >= ZERO && N*TOL*( SMIN / SMAX ) <= MAX( EPS, HNDRTH*TOL ) ) {

         // Use a zero shift to avoid loss of relative accuracy

         SHIFT = ZERO;
      } else {

         // Compute the shift from 2-by-2 block at end of matrix

         if ( IDIR == 1 ) {
            SLL = ABS( D( LL ) );
            slas2(D( M-1 ), E( M-1 ), D( M ), SHIFT, R );
         } else {
            SLL = ABS( D( M ) );
            slas2(D( LL ), E( LL ), D( LL+1 ), SHIFT, R );
         }

         // Test if shift negligible, and if so set to zero

         if ( SLL > ZERO ) {
            IF( ( SHIFT / SLL )**2 < EPS ) SHIFT = ZERO;
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
               slartg(D( I )*CS, E( I ), CS, SN, R );
               if (I > LL) E( I-1 ) = OLDSN*R;
               slartg(OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) );
               RWORK( I-LL+1 ) = CS;
               RWORK( I-LL+1+NM1 ) = SN;
               RWORK( I-LL+1+NM12 ) = OLDCS;
               RWORK( I-LL+1+NM13 ) = OLDSN;
            } // 120
            H = D( M )*CS;
            D( M ) = H*OLDCS;
            E( M-1 ) = H*OLDSN;

            // Update singular vectors

            if (NCVT > 0) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ), RWORK( N ), VT( LL, 1 ), LDVT )             IF( NRU > 0 ) CALL CLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ), RWORK( NM13+1 ), U( 1, LL ), LDU )             IF( NCC > 0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ), RWORK( NM13+1 ), C( LL, 1 ), LDC );

            // Test convergence

            IF( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = ZERO;

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            CS = ONE;
            OLDCS = ONE;
            DO 130 I = M, LL + 1, -1;
               slartg(D( I )*CS, E( I-1 ), CS, SN, R );
               if (I < M) E( I ) = OLDSN*R;
               slartg(OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) );
               RWORK( I-LL ) = CS;
               RWORK( I-LL+NM1 ) = -SN;
               RWORK( I-LL+NM12 ) = OLDCS;
               RWORK( I-LL+NM13 ) = -OLDSN;
            } // 130
            H = D( LL )*CS;
            D( LL ) = H*OLDCS;
            E( LL ) = H*OLDSN;

            // Update singular vectors

            if (NCVT > 0) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ), RWORK( NM13+1 ), VT( LL, 1 ), LDVT )             IF( NRU > 0 ) CALL CLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ), RWORK( N ), U( 1, LL ), LDU )             IF( NCC > 0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ), RWORK( N ), C( LL, 1 ), LDC );

            // Test convergence

            IF( ABS( E( LL ) ) <= THRESH ) E( LL ) = ZERO;
         }
      } else {

         // Use nonzero shift

         if ( IDIR == 1 ) {

            // Chase bulge from top to bottom
            // Save cosines and sines for later singular vector updates

            F = ( ABS( D( LL ) )-SHIFT )* ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) );
            G = E( LL );
            for (I = LL; I <= M - 1; I++) { // 140
               slartg(F, G, COSR, SINR, R );
               if (I > LL) E( I-1 ) = R;
               F = COSR*D( I ) + SINR*E( I );
               E( I ) = COSR*E( I ) - SINR*D( I );
               G = SINR*D( I+1 );
               D( I+1 ) = COSR*D( I+1 );
               slartg(F, G, COSL, SINL, R );
               D( I ) = R;
               F = COSL*E( I ) + SINL*D( I+1 );
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I );
               if ( I < M-1 ) {
                  G = SINL*E( I+1 );
                  E( I+1 ) = COSL*E( I+1 );
               }
               RWORK( I-LL+1 ) = COSR;
               RWORK( I-LL+1+NM1 ) = SINR;
               RWORK( I-LL+1+NM12 ) = COSL;
               RWORK( I-LL+1+NM13 ) = SINL;
            } // 140
            E( M-1 ) = F;

            // Update singular vectors

            if (NCVT > 0) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ), RWORK( N ), VT( LL, 1 ), LDVT )             IF( NRU > 0 ) CALL CLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ), RWORK( NM13+1 ), U( 1, LL ), LDU )             IF( NCC > 0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ), RWORK( NM13+1 ), C( LL, 1 ), LDC );

            // Test convergence

            IF( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = ZERO;

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / D( M ) );
            G = E( M-1 );
            DO 150 I = M, LL + 1, -1;
               slartg(F, G, COSR, SINR, R );
               if (I < M) E( I ) = R;
               F = COSR*D( I ) + SINR*E( I-1 );
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I );
               G = SINR*D( I-1 );
               D( I-1 ) = COSR*D( I-1 );
               slartg(F, G, COSL, SINL, R );
               D( I ) = R;
               F = COSL*E( I-1 ) + SINL*D( I-1 );
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 );
               if ( I > LL+1 ) {
                  G = SINL*E( I-2 );
                  E( I-2 ) = COSL*E( I-2 );
               }
               RWORK( I-LL ) = COSR;
               RWORK( I-LL+NM1 ) = -SINR;
               RWORK( I-LL+NM12 ) = COSL;
               RWORK( I-LL+NM13 ) = -SINL;
            } // 150
            E( LL ) = F;

            // Test convergence

            IF( ABS( E( LL ) ) <= THRESH ) E( LL ) = ZERO;

            // Update singular vectors if desired

            if (NCVT > 0) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ), RWORK( NM13+1 ), VT( LL, 1 ), LDVT )             IF( NRU > 0 ) CALL CLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ), RWORK( N ), U( 1, LL ), LDU )             IF( NCC > 0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ), RWORK( N ), C( LL, 1 ), LDC );
         }
      }

      // QR iteration finished, go back and check convergence

      GO TO 60;

      // All singular values converged, so make them positive

      } // 160
      for (I = 1; I <= N; I++) { // 170
         if ( D( I ) < ZERO ) {
            D( I ) = -D( I );

            // Change sign of singular vectors, if desired

            if (NCVT > 0) CALL CSSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT );
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

            D( ISUB ) = D( N+1-I );
            D( N+1-I ) = SMIN;
            if (NCVT > 0) CALL CSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), LDVT );
            if (NRU > 0) CALL CSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )             IF( NCC > 0 ) CALL CSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC );
         }
      } // 190
      GO TO 220;

      // Maximum number of iterations exceeded, failure to converge

      } // 200
      INFO = 0;
      for (I = 1; I <= N - 1; I++) { // 210
         IF( E( I ) != ZERO ) INFO = INFO + 1;
      } // 210
      } // 220
      return;

      // End of CBDSQR

      }
