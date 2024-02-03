      SUBROUTINE CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            C( LDC, * ), U( LDU, * ), VT( LDVT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               NEGONE
      const              NEGONE = -1.0E0 ;
      REAL               HNDRTH
      const              HNDRTH = 0.01E0 ;
      REAL               TEN
      const              TEN = 10.0E0 ;
      REAL               HNDRD
      const              HNDRD = 100.0E0 ;
      REAL               MEIGTH
      const              MEIGTH = -0.125E0 ;
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
      REAL               SLAMCH
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

      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NCVT.LT.0 ) {
         INFO = -3
      } else if ( NRU.LT.0 ) {
         INFO = -4
      } else if ( NCC.LT.0 ) {
         INFO = -5
      } else if ( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR. ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) {
         INFO = -9
      } else if ( LDU.LT.MAX( 1, NRU ) ) {
         INFO = -11
      } else if ( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR. ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) {
         INFO = -13
      }
      if ( INFO.NE.0 ) {
         xerbla('CBDSQR', -INFO );
         RETURN
      }
      IF( N.EQ.0 ) RETURN       IF( N.EQ.1 ) GO TO 160

      // ROTATE is true if any singular vectors desired, false otherwise

      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )

      // If no singular vectors desired, use qd algorithm

      if ( .NOT.ROTATE ) {
         slasq1(N, D, E, RWORK, INFO );

      // If INFO equals 2, dqds didn't finish, try to finish

         IF( INFO .NE. 2 ) RETURN
         INFO = 0
      }

      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0

      // Get machine constants

      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      if ( LOWER ) {
         DO 10 I = 1, N - 1
            slartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            RWORK( I ) = CS
            RWORK( NM1+I ) = SN
   10    CONTINUE

         // Update singular vectors if desired

         IF( NRU.GT.0 ) CALL CLASR( 'R', 'V', 'F', NRU, N, RWORK( 1 ), RWORK( N ), U, LDU )          IF( NCC.GT.0 ) CALL CLASR( 'L', 'V', 'F', N, NCC, RWORK( 1 ), RWORK( N ), C, LDC )
      }

      // Compute singular values to relative accuracy TOL
      // (By setting TOL to be negative, algorithm will compute
      // singular values to absolute accuracy ABS(TOL)*norm(input matrix))

      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS

      // Compute approximate maximum, minimum singular values

      SMAX = ZERO
      DO 20 I = 1, N
         SMAX = MAX( SMAX, ABS( D( I ) ) )
   20 CONTINUE
      DO 30 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( E( I ) ) )
   30 CONTINUE
      SMIN = ZERO
      if ( TOL.GE.ZERO ) {

         // Relative accuracy desired

         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO ) GO TO 50
         MU = SMINOA
         DO 40 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO ) GO TO 50
   40    CONTINUE
   50    CONTINUE
         SMINOA = SMINOA / SQRT( REAL( N ) )
         THRESH = MAX( TOL*SMINOA, MAXITR*(N*(N*UNFL)) )
      } else {

         // Absolute accuracy desired

         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*(N*(N*UNFL)) )
      }

      // Prepare for main iteration loop for the singular values
      // (MAXIT is the maximum number of passes through the inner
      // loop permitted before nonconvergence signalled.)

      MAXITDIVN = MAXITR*N
      ITERDIVN = 0
      ITER = -1
      OLDLL = -1
      OLDM = -1

      // M points to last element of unconverged part of matrix

      M = N

      // Begin main iteration loop

   60 CONTINUE

      // Check for convergence or exceeding iteration count

      IF( M.LE.1 ) GO TO 160
      if ( ITER.GE.N ) {
         ITER = ITER - N
         ITERDIVN = ITERDIVN + 1
         IF( ITERDIVN.GE.MAXITDIVN ) GO TO 200
      }

      // Find diagonal block of matrix to work on

      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH ) D( M ) = ZERO
      SMAX = ABS( D( M ) )
      DO 70 LLL = 1, M - 1
         LL = M - LLL
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH ) D( LL ) = ZERO          IF( ABSE.LE.THRESH ) GO TO 80
         SMAX = MAX( SMAX, ABSS, ABSE )
   70 CONTINUE
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO

      // Matrix splits since E(LL) = 0

      if ( LL.EQ.M-1 ) {

         // Convergence of bottom singular value, return to top of loop

         M = M - 1
         GO TO 60
      }
   90 CONTINUE
      LL = LL + 1

      // E(LL) through E(M-1) are nonzero, E(LL-1) is zero

      if ( LL.EQ.M-1 ) {

         // 2 by 2 block, handle separately

         slasv2(D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, COSR, SINL, COSL );
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN

         // Compute singular vectors, if desired

         IF( NCVT.GT.0 ) CALL CSROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, SINR )          IF( NRU.GT.0 ) CALL CSROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )          IF( NCC.GT.0 ) CALL CSROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, SINL )
         M = M - 2
         GO TO 60
      }

      // If working on new submatrix, choose shift direction
      // (from larger end diagonal element towards smaller)

      if ( LL.GT.OLDM .OR. M.LT.OLDLL ) {
         if ( ABS( D( LL ) ).GE.ABS( D( M ) ) ) {

            // Chase bulge from top (big end) to bottom (small end)

            IDIR = 1
         } else {

            // Chase bulge from bottom (big end) to top (small end)

            IDIR = 2
         }
      }

      // Apply convergence tests

      if ( IDIR.EQ.1 ) {

         // Run convergence test in forward direction
         // First apply standard test to bottom of matrix

         if ( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR. ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) {
            E( M-1 ) = ZERO
            GO TO 60
         }

         if ( TOL.GE.ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion forward

            MU = ABS( D( LL ) )
            SMIN = MU
            DO 100 LLL = LL, M - 1
               if ( ABS( E( LLL ) ).LE.TOL*MU ) {
                  E( LLL ) = ZERO
                  GO TO 60
               }
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMIN = MIN( SMIN, MU )
  100       CONTINUE
         }

      } else {

         // Run convergence test in backward direction
         // First apply standard test to top of matrix

         if ( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR. ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) {
            E( LL ) = ZERO
            GO TO 60
         }

         if ( TOL.GE.ZERO ) {

            // If relative accuracy desired,
            // apply convergence criterion backward

            MU = ABS( D( M ) )
            SMIN = MU
            DO 110 LLL = M - 1, LL, -1
               if ( ABS( E( LLL ) ).LE.TOL*MU ) {
                  E( LLL ) = ZERO
                  GO TO 60
               }
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMIN = MIN( SMIN, MU )
  110       CONTINUE
         }
      }
      OLDLL = LL
      OLDM = M

      // Compute shift.  First, test if shifting would ruin relative
      // accuracy, and if so set the shift to zero.

      if ( TOL.GE.ZERO .AND. N*TOL*( SMIN / SMAX ).LE. MAX( EPS, HNDRTH*TOL ) ) {

         // Use a zero shift to avoid loss of relative accuracy

         SHIFT = ZERO
      } else {

         // Compute the shift from 2-by-2 block at end of matrix

         if ( IDIR.EQ.1 ) {
            SLL = ABS( D( LL ) )
            slas2(D( M-1 ), E( M-1 ), D( M ), SHIFT, R );
         } else {
            SLL = ABS( D( M ) )
            slas2(D( LL ), E( LL ), D( LL+1 ), SHIFT, R );
         }

         // Test if shift negligible, and if so set to zero

         if ( SLL.GT.ZERO ) {
            IF( ( SHIFT / SLL )**2.LT.EPS ) SHIFT = ZERO
         }
      }

      // Increment iteration count

      ITER = ITER + M - LL

      // If SHIFT = 0, do simplified QR iteration

      if ( SHIFT.EQ.ZERO ) {
         if ( IDIR.EQ.1 ) {

            // Chase bulge from top to bottom
            // Save cosines and sines for later singular vector updates

            CS = ONE
            OLDCS = ONE
            DO 120 I = LL, M - 1
               slartg(D( I )*CS, E( I ), CS, SN, R );
               IF( I.GT.LL ) E( I-1 ) = OLDSN*R
               slartg(OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) );
               RWORK( I-LL+1 ) = CS
               RWORK( I-LL+1+NM1 ) = SN
               RWORK( I-LL+1+NM12 ) = OLDCS
               RWORK( I-LL+1+NM13 ) = OLDSN
  120       CONTINUE
            H = D( M )*CS
            D( M ) = H*OLDCS
            E( M-1 ) = H*OLDSN

            // Update singular vectors

            IF( NCVT.GT.0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ), RWORK( N ), VT( LL, 1 ), LDVT )             IF( NRU.GT.0 ) CALL CLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ), RWORK( NM13+1 ), U( 1, LL ), LDU )             IF( NCC.GT.0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ), RWORK( NM13+1 ), C( LL, 1 ), LDC )

            // Test convergence

            IF( ABS( E( M-1 ) ).LE.THRESH ) E( M-1 ) = ZERO

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            CS = ONE
            OLDCS = ONE
            DO 130 I = M, LL + 1, -1
               slartg(D( I )*CS, E( I-1 ), CS, SN, R );
               IF( I.LT.M ) E( I ) = OLDSN*R
               slartg(OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) );
               RWORK( I-LL ) = CS
               RWORK( I-LL+NM1 ) = -SN
               RWORK( I-LL+NM12 ) = OLDCS
               RWORK( I-LL+NM13 ) = -OLDSN
  130       CONTINUE
            H = D( LL )*CS
            D( LL ) = H*OLDCS
            E( LL ) = H*OLDSN

            // Update singular vectors

            IF( NCVT.GT.0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ), RWORK( NM13+1 ), VT( LL, 1 ), LDVT )             IF( NRU.GT.0 ) CALL CLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ), RWORK( N ), U( 1, LL ), LDU )             IF( NCC.GT.0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ), RWORK( N ), C( LL, 1 ), LDC )

            // Test convergence

            IF( ABS( E( LL ) ).LE.THRESH ) E( LL ) = ZERO
         }
      } else {

         // Use nonzero shift

         if ( IDIR.EQ.1 ) {

            // Chase bulge from top to bottom
            // Save cosines and sines for later singular vector updates

            F = ( ABS( D( LL ) )-SHIFT )* ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
            DO 140 I = LL, M - 1
               slartg(F, G, COSR, SINR, R );
               IF( I.GT.LL ) E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               slartg(F, G, COSL, SINL, R );
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               if ( I.LT.M-1 ) {
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
               }
               RWORK( I-LL+1 ) = COSR
               RWORK( I-LL+1+NM1 ) = SINR
               RWORK( I-LL+1+NM12 ) = COSL
               RWORK( I-LL+1+NM13 ) = SINL
  140       CONTINUE
            E( M-1 ) = F

            // Update singular vectors

            IF( NCVT.GT.0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ), RWORK( N ), VT( LL, 1 ), LDVT )             IF( NRU.GT.0 ) CALL CLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ), RWORK( NM13+1 ), U( 1, LL ), LDU )             IF( NCC.GT.0 ) CALL CLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ), RWORK( NM13+1 ), C( LL, 1 ), LDC )

            // Test convergence

            IF( ABS( E( M-1 ) ).LE.THRESH ) E( M-1 ) = ZERO

         } else {

            // Chase bulge from bottom to top
            // Save cosines and sines for later singular vector updates

            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / D( M ) )
            G = E( M-1 )
            DO 150 I = M, LL + 1, -1
               slartg(F, G, COSR, SINR, R );
               IF( I.LT.M ) E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               slartg(F, G, COSL, SINL, R );
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               if ( I.GT.LL+1 ) {
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
               }
               RWORK( I-LL ) = COSR
               RWORK( I-LL+NM1 ) = -SINR
               RWORK( I-LL+NM12 ) = COSL
               RWORK( I-LL+NM13 ) = -SINL
  150       CONTINUE
            E( LL ) = F

            // Test convergence

            IF( ABS( E( LL ) ).LE.THRESH ) E( LL ) = ZERO

            // Update singular vectors if desired

            IF( NCVT.GT.0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ), RWORK( NM13+1 ), VT( LL, 1 ), LDVT )             IF( NRU.GT.0 ) CALL CLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ), RWORK( N ), U( 1, LL ), LDU )             IF( NCC.GT.0 ) CALL CLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ), RWORK( N ), C( LL, 1 ), LDC )
         }
      }

      // QR iteration finished, go back and check convergence

      GO TO 60

      // All singular values converged, so make them positive

  160 CONTINUE
      DO 170 I = 1, N
         if ( D( I ).LT.ZERO ) {
            D( I ) = -D( I )

            // Change sign of singular vectors, if desired

            IF( NCVT.GT.0 ) CALL CSSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         }
  170 CONTINUE

      // Sort the singular values into decreasing order (insertion sort on
      // singular values, but only one transposition per singular vector)

      DO 190 I = 1, N - 1

         // Scan for smallest D(I)

         ISUB = 1
         SMIN = D( 1 )
         DO 180 J = 2, N + 1 - I
            if ( D( J ).LE.SMIN ) {
               ISUB = J
               SMIN = D( J )
            }
  180    CONTINUE
         if ( ISUB.NE.N+1-I ) {

            // Swap singular values and vectors

            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 ) CALL CSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), LDVT )
            IF( NRU.GT.0 ) CALL CSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )             IF( NCC.GT.0 ) CALL CSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         }
  190 CONTINUE
      GO TO 220

      // Maximum number of iterations exceeded, failure to converge

  200 CONTINUE
      INFO = 0
      DO 210 I = 1, N - 1
         IF( E( I ).NE.ZERO ) INFO = INFO + 1
  210 CONTINUE
  220 CONTINUE
      RETURN

      // End of CBDSQR

      }
