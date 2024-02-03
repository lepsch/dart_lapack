      void dgeevx(BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, INFO ) {
      // implicit none

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N;
      double             ABNRM;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), RCONDE( * ), RCONDV( * ), SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WORK( * ), WR( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
      String             JOB, SIDE;
      int                HSWORK, I, ICOND, IERR, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double             ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, SN;
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC3, DTRSNA, XERBLA
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                IDAMAX, ILAENV;
      //- double             DLAMCH, DLANGE, DLAPY2, DNRM2;
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
      WNTSNN = LSAME( SENSE, 'N' );
      WNTSNE = LSAME( SENSE, 'E' );
      WNTSNV = LSAME( SENSE, 'V' );
      WNTSNB = LSAME( SENSE, 'B' );
      if ( !( LSAME( BALANC, 'N' ) || LSAME( BALANC, 'S' ) || LSAME( BALANC, 'P' ) || LSAME( BALANC, 'B' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTVL ) && ( !LSAME( JOBVL, 'N' ) ) ) {
         INFO = -2;
      } else if ( ( !WANTVR ) && ( !LSAME( JOBVR, 'N' ) ) ) {
         INFO = -3;
      } else if ( !( WNTSNN || WNTSNE || WNTSNB || WNTSNV ) || ( ( WNTSNE || WNTSNB ) && !( WANTVL && WANTVR ) ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDVL < 1 || ( WANTVL && LDVL < N ) ) {
         INFO = -11;
      } else if ( LDVR < 1 || ( WANTVR && LDVR < N ) ) {
         INFO = -13;
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
            MAXWRK = N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 );

            if ( WANTVL ) {
               dtrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO );
            } else if ( WANTVR ) {
               dtrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               dhseqr('S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO );
            } else {
               if ( WNTSNN ) {
                  dhseqr('E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO );
               } else {
                  dhseqr('S', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO );
               }
            }
            HSWORK = INT( WORK(1) );

            if ( ( !WANTVL ) && ( !WANTVR ) ) {
               MINWRK = 2*N;
               if ( !WNTSNN) MINWRK = max( MINWRK, N*N+6*N );
               MAXWRK = max( MAXWRK, HSWORK );
               if ( !WNTSNN) MAXWRK = max( MAXWRK, N*N + 6*N );
            } else {
               MINWRK = 3*N;
               if( ( !WNTSNN ) && ( !WNTSNE ) ) MINWRK = max( MINWRK, N*N + 6*N );
               MAXWRK = max( MAXWRK, HSWORK );
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )                IF( ( !WNTSNN ) && ( !WNTSNE ) ) MAXWRK = max( MAXWRK, N*N + 6*N );
               MAXWRK = max( MAXWRK, 3*N );
            }
            MAXWRK = max( MAXWRK, MINWRK );
         }
         WORK( 1 ) = MAXWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -21;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGEEVX', -INFO );
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

      ICOND = 0;
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM );
      SCALEA = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         SCALEA = true;
         CSCALE = SMLNUM;
      } else if ( ANRM > BIGNUM ) {
         SCALEA = true;
         CSCALE = BIGNUM;
      }
      if (SCALEA) dlascl( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR );

      // Balance the matrix and compute ABNRM

      dgebal(BALANC, N, A, LDA, ILO, IHI, SCALE, IERR );
      ABNRM = DLANGE( '1', N, N, A, LDA, DUM );
      if ( SCALEA ) {
         DUM( 1 ) = ABNRM;
         dlascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR );
         ABNRM = DUM( 1 );
      }

      // Reduce to upper Hessenberg form
      // (Workspace: need 2*N, prefer N+N*NB)

      ITAU = 1;
      IWRK = ITAU + N;
      dgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L';
         dlacpy('L', N, N, A, LDA, VL, LDVL );

         // Generate orthogonal matrix in VL
         // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

         dorghr(N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VL
         // (Workspace: need 1, prefer HSWORK (see comments) )

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
         // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

         dorghr(N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VR
         // (Workspace: need 1, prefer HSWORK (see comments) )

         IWRK = ITAU;
         dhseqr('S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );

      } else {

         // Compute eigenvalues only
         // If condition numbers desired, compute Schur form

         if ( WNTSNN ) {
            JOB = 'E';
         } else {
            JOB = 'S';
         }

         // (Workspace: need 1, prefer HSWORK (see comments) )

         IWRK = ITAU;
         dhseqr(JOB, 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );
      }

      // If INFO != 0 from DHSEQR, then quit

      if (INFO != 0) GO TO 50;

      if ( WANTVL || WANTVR ) {

         // Compute left and/or right eigenvectors
         // (Workspace: need 3*N, prefer N + 2*N*NB)

         dtrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR );
      }

      // Compute condition numbers if desired
      // (Workspace: need N*N+6*N unless SENSE = 'E')

      if ( !WNTSNN ) {
         dtrsna(SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, IWORK, ICOND );
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors

         dgebak(BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, IERR );

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
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2;
               } // 10
               K = IDAMAX( N, WORK, 1 );
               dlartg(VL( K, I ), VL( K, I+1 ), CS, SN, R );
               drot(N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN );
               VL( K, I+1 ) = ZERO;
            }
         } // 20
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors

         dgebak(BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, IERR );

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
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2;
               } // 30
               K = IDAMAX( N, WORK, 1 );
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
         if ( INFO == 0 ) {
            if( ( WNTSNV || WNTSNB ) && ICOND == 0 ) dlascl( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, IERR );
         } else {
            dlascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, IERR );
            dlascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, IERR );
         }
      }

      WORK( 1 ) = MAXWRK;
      return;
      }
