      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO );
      // implicit none

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WORK( * ), WR( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR;
      String             SIDE;
      int                HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double             ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, SN;
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC3, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX, ILAENV;
      double             DLAMCH, DLANGE, DLAPY2, DNRM2;
      // EXTERNAL LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      WANTVL = LSAME( JOBVL, 'V' );
      WANTVR = LSAME( JOBVR, 'V' );
      if ( ( !WANTVL ) && ( !LSAME( JOBVL, 'N' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTVR ) && ( !LSAME( JOBVR, 'N' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDVL < 1 || ( WANTVL && LDVL < N ) ) {
         INFO = -9;
      } else if ( LDVR < 1 || ( WANTVR && LDVR < N ) ) {
         INFO = -11;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by DHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            MAXWRK = 1;
         } else {
            MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 );
            if ( WANTVL ) {
               MINWRK = 4*N;
               MAXWRK = max( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) );
               dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO );
               HSWORK = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + 1, N + HSWORK );
               dtrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               MAXWRK = max( MAXWRK, 4*N );
            } else if ( WANTVR ) {
               MINWRK = 4*N;
               MAXWRK = max( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) );
               dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO );
               HSWORK = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + 1, N + HSWORK );
               dtrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               MAXWRK = max( MAXWRK, 4*N );
            } else {
               MINWRK = 3*N;
               dhseqr('E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO );
               HSWORK = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + 1, N + HSWORK );
            }
            MAXWRK = max( MAXWRK, MINWRK );
         }
         WORK( 1 ) = MAXWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -13;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGEEV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, DUM );
      SCALEA = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         SCALEA = true;
         CSCALE = SMLNUM;
      } else if ( ANRM > BIGNUM ) {
         SCALEA = true;
         CSCALE = BIGNUM;
      }
      if (SCALEA) CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR );

      // Balance the matrix
      // (Workspace: need N)

      IBAL = 1;
      dgebal('B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR );

      // Reduce to upper Hessenberg form
      // (Workspace: need 3*N, prefer 2*N+N*NB)

      ITAU = IBAL + N;
      IWRK = ITAU + N;
      dgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L';
         dlacpy('L', N, N, A, LDA, VL, LDVL );

         // Generate orthogonal matrix in VL
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         dorghr(N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VL
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU;
         dhseqr('S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO );

         if ( WANTVR ) {

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B';
            dlacpy('F', N, N, VL, LDVL, VR, LDVR );
         }

      } else if ( WANTVR ) {

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R';
         dlacpy('L', N, N, A, LDA, VR, LDVR );

         // Generate orthogonal matrix in VR
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         dorghr(N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VR
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU;
         dhseqr('S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );

      } else {

         // Compute eigenvalues only
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU;
         dhseqr('E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );
      }

      // If INFO != 0 from DHSEQR, then quit

      if (INFO != 0) GO TO 50;

      if ( WANTVL || WANTVR ) {

         // Compute left and/or right eigenvectors
         // (Workspace: need 4*N, prefer N + N + 2*N*NB)

         dtrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR );
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors
         // (Workspace: need N)

         dgebak('B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, IERR );

         // Normalize left eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 20
            if ( WI( I ) == ZERO ) {
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 );
               dscal(N, SCL, VL( 1, I ), 1 );
            } else if ( WI( I ) > ZERO ) {
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ), DNRM2( N, VL( 1, I+1 ), 1 ) );
               dscal(N, SCL, VL( 1, I ), 1 );
               dscal(N, SCL, VL( 1, I+1 ), 1 );
               for (K = 1; K <= N; K++) { // 10
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2;
               } // 10
               K = IDAMAX( N, WORK( IWRK ), 1 );
               dlartg(VL( K, I ), VL( K, I+1 ), CS, SN, R );
               drot(N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN );
               VL( K, I+1 ) = ZERO;
            }
         } // 20
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors
         // (Workspace: need N)

         dgebak('B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, IERR );

         // Normalize right eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 40
            if ( WI( I ) == ZERO ) {
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 );
               dscal(N, SCL, VR( 1, I ), 1 );
            } else if ( WI( I ) > ZERO ) {
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ), DNRM2( N, VR( 1, I+1 ), 1 ) );
               dscal(N, SCL, VR( 1, I ), 1 );
               dscal(N, SCL, VR( 1, I+1 ), 1 );
               for (K = 1; K <= N; K++) { // 30
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2;
               } // 30
               K = IDAMAX( N, WORK( IWRK ), 1 );
               dlartg(VR( K, I ), VR( K, I+1 ), CS, SN, R );
               drot(N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN );
               VR( K, I+1 ) = ZERO;
            }
         } // 40
      }

      // Undo scaling if necessary

      } // 50
      if ( SCALEA ) {
         dlascl('G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), max( N-INFO, 1 ), IERR );
         dlascl('G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), max( N-INFO, 1 ), IERR );
         if ( INFO > 0 ) {
            dlascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, IERR );
            dlascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, IERR );
         }
      }

      WORK( 1 ) = MAXWRK;
      return;

      // End of DGEEV

      }
