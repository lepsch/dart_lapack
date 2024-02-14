      void cgeevx(final int BALANC, final int JOBVL, final int JOBVR, final int SENSE, final int N, final Matrix<double> A_, final int LDA, final int W, final Matrix<double> VL_, final int LDVL, final Matrix<double> VR_, final int LDVR, final int ILO, final int IHI, final int SCALE, final int ABNRM, final int RCONDE, final int RCONDV, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Box<int> INFO,) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim();
  final VL = VL_.dim();
  final VR = VR_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N;
      double               ABNRM;
      double               RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )       Complex            A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
      String             JOB, SIDE;
      int                HSWORK, I, ICOND, IERR, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double               ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM;
      Complex            TMP;
      bool               SELECT( 1 );
      double   DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASCL, XERBLA, CSSCAL, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY, CLASCL, CSCAL, CTREVC3, CTRSNA, CUNGHR
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX, ILAENV;
      //- REAL               SLAMCH, SCNRM2, CLANGE, SROUNDUP_LWORK;
      // EXTERNAL lsame, ISAMAX, ILAENV, SLAMCH, SCNRM2, CLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, CMPLX, CONJG, AIMAG, MAX, SQRT

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      WANTVL = lsame( JOBVL, 'V' );
      WANTVR = lsame( JOBVR, 'V' );
      WNTSNN = lsame( SENSE, 'N' );
      WNTSNE = lsame( SENSE, 'E' );
      WNTSNV = lsame( SENSE, 'V' );
      WNTSNB = lsame( SENSE, 'B' );
      if ( !( lsame( BALANC, 'N' ) || lsame( BALANC, 'S' ) || lsame( BALANC, 'P' ) || lsame( BALANC, 'B' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTVL ) && ( !lsame( JOBVL, 'N' ) ) ) {
         INFO = -2;
      } else if ( ( !WANTVR ) && ( !lsame( JOBVR, 'N' ) ) ) {
         INFO = -3;
      } else if ( !( WNTSNN || WNTSNE || WNTSNB || WNTSNV ) || ( ( WNTSNE || WNTSNB ) && !( WANTVL && WANTVR ) ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDVL < 1 || ( WANTVL && LDVL < N ) ) {
         INFO = -10;
      } else if ( LDVR < 1 || ( WANTVR && LDVR < N ) ) {
         INFO = -12;
      }

      // Compute workspace
      //  (Note: Comments in the code beginning "Workspace:" describe the
      //   minimal amount of workspace needed at that point in the code,
      //   as well as the preferred amount for good performance.
      //   CWorkspace refers to complex workspace, and RWorkspace to real
      //   workspace. NB refers to the optimal block size for the
      //   immediately following subroutine, as returned by ILAENV.
      //   HSWORK refers to the workspace preferred by CHSEQR, as
      //   calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
      //   the worst case.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            MAXWRK = 1;
         } else {
            MAXWRK = N + N*ilaenv( 1, 'CGEHRD', ' ', N, 1, N, 0 );

            if ( WANTVL ) {
               ctrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, LWORK_TREVC );
               chseqr('S', 'V', N, 1, N, A, LDA, W, VL, LDVL, WORK, -1, INFO );
            } else if ( WANTVR ) {
               ctrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, LWORK_TREVC );
               chseqr('S', 'V', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO );
            } else {
               if ( WNTSNN ) {
                  chseqr('E', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO );
               } else {
                  chseqr('S', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO );
               }
            }
            HSWORK = INT( WORK(1) );

            if ( ( !WANTVL ) && ( !WANTVR ) ) {
               MINWRK = 2*N;
               if( !( WNTSNN || WNTSNE ) ) MINWRK = max( MINWRK, N*N + 2*N );
               MAXWRK = max( MAXWRK, HSWORK );
               if( !( WNTSNN || WNTSNE ) ) MAXWRK = max( MAXWRK, N*N + 2*N );
            } else {
               MINWRK = 2*N;
               if( !( WNTSNN || WNTSNE ) ) MINWRK = max( MINWRK, N*N + 2*N );
               MAXWRK = max( MAXWRK, HSWORK );
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ilaenv( 1, 'CUNGHR', ' ', N, 1, N, -1 ) )                IF( !( WNTSNN || WNTSNE ) ) MAXWRK = max( MAXWRK, N*N + 2*N );
               MAXWRK = max( MAXWRK, 2*N );
            }
            MAXWRK = max( MAXWRK, MINWRK );
         }
         WORK[1] = SROUNDUP_LWORK(MAXWRK);

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -20;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGEEVX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ICOND = 0;
      ANRM = CLANGE( 'M', N, N, A, LDA, DUM );
      SCALEA = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         SCALEA = true;
         CSCALE = SMLNUM;
      } else if ( ANRM > BIGNUM ) {
         SCALEA = true;
         CSCALE = BIGNUM;
      }
      if (SCALEA) clascl( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR );

      // Balance the matrix and compute ABNRM

      cgebal(BALANC, N, A, LDA, ILO, IHI, SCALE, IERR );
      ABNRM = CLANGE( '1', N, N, A, LDA, DUM );
      if ( SCALEA ) {
         DUM[1] = ABNRM;
         slascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR );
         ABNRM = DUM( 1 );
      }

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1;
      IWRK = ITAU + N;
      cgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L';
         clacpy('L', N, N, A, LDA, VL, LDVL );

         // Generate unitary matrix in VL
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         cunghr(N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VL
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         chseqr('S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO );

         if ( WANTVR ) {

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B';
            clacpy('F', N, N, VL, LDVL, VR, LDVR );
         }

      } else if ( WANTVR ) {

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R';
         clacpy('L', N, N, A, LDA, VR, LDVR );

         // Generate unitary matrix in VR
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         cunghr(N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VR
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         chseqr('S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );

      } else {

         // Compute eigenvalues only
         // If condition numbers desired, compute Schur form

         if ( WNTSNN ) {
            JOB = 'E';
         } else {
            JOB = 'S';
         }

         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         chseqr(JOB, 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );
      }

      // If INFO != 0 from CHSEQR, then quit

      if (INFO != 0) GO TO 50;

      if ( WANTVL || WANTVR ) {

         // Compute left and/or right eigenvectors
         // (CWorkspace: need 2*N, prefer N + 2*N*NB)
         // (RWorkspace: need N)

         ctrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, RWORK, N, IERR );
      }

      // Compute condition numbers if desired
      // (CWorkspace: need N*N+2*N unless SENSE = 'E')
      // (RWorkspace: need 2*N unless SENSE = 'E')

      if ( !WNTSNN ) {
         ctrsna(SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, RWORK, ICOND );
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors

         cgebak(BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, IERR );

         // Normalize left eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 20
            SCL = ONE / SCNRM2( N, VL( 1, I ), 1 );
            csscal(N, SCL, VL( 1, I ), 1 );
            for (K = 1; K <= N; K++) { // 10
               RWORK[K] = double( VL( K, I ) )**2 + AIMAG( VL( K, I ) )**2;
            } // 10
            K = ISAMAX( N, RWORK, 1 );
            TMP = CONJG( VL( K, I ) ) / sqrt( RWORK( K ) );
            cscal(N, TMP, VL( 1, I ), 1 );
            VL[K][I] = CMPLX( double( VL( K, I ) ), ZERO );
         } // 20
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors

         cgebak(BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, IERR );

         // Normalize right eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 40
            SCL = ONE / SCNRM2( N, VR( 1, I ), 1 );
            csscal(N, SCL, VR( 1, I ), 1 );
            for (K = 1; K <= N; K++) { // 30
               RWORK[K] = double( VR( K, I ) )**2 + AIMAG( VR( K, I ) )**2;
            } // 30
            K = ISAMAX( N, RWORK, 1 );
            TMP = CONJG( VR( K, I ) ) / sqrt( RWORK( K ) );
            cscal(N, TMP, VR( 1, I ), 1 );
            VR[K][I] = CMPLX( double( VR( K, I ) ), ZERO );
         } // 40
      }

      // Undo scaling if necessary

      } // 50
      if ( SCALEA ) {
         clascl('G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), max( N-INFO, 1 ), IERR );
         if ( INFO == 0 ) {
            if( ( WNTSNV || WNTSNB ) && ICOND == 0 ) slascl( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, IERR );
         } else {
            clascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR );
         }
      }

      WORK[1] = SROUNDUP_LWORK(MAXWRK);
      }
