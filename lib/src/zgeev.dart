      void zgeev(final int JOBVL, final int JOBVR, final int N, final Matrix<double> A, final int LDA, final int W, final Matrix<double> VL, final int LDVL, final Matrix<double> VR, final int LDVR, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Box<int> INFO,) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDVL, LDVR, LWORK, N;
      double             RWORK( * );
      Complex         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, SCALEA, WANTVL, WANTVR;
      String             SIDE;
      int                HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double             ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM;
      Complex         TMP;
      bool               SELECT( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC3, ZUNGHR
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax, ILAENV;
      //- double             DLAMCH, DZNRM2, ZLANGE;
      // EXTERNAL lsame, idamax, ILAENV, DLAMCH, DZNRM2, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, CONJG, AIMAG, MAX, SQRT

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      WANTVL = lsame( JOBVL, 'V' );
      WANTVR = lsame( JOBVR, 'V' );
      if ( ( !WANTVL ) && ( !lsame( JOBVL, 'N' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTVR ) && ( !lsame( JOBVR, 'N' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDVL < 1 || ( WANTVL && LDVL < N ) ) {
         INFO = -8;
      } else if ( LDVR < 1 || ( WANTVR && LDVR < N ) ) {
         INFO = -10;
      }

      // Compute workspace
      //  (Note: Comments in the code beginning "Workspace:" describe the
      //   minimal amount of workspace needed at that point in the code,
      //   as well as the preferred amount for good performance.
      //   CWorkspace refers to complex workspace, and RWorkspace to real
      //   workspace. NB refers to the optimal block size for the
      //   immediately following subroutine, as returned by ILAENV.
      //   HSWORK refers to the workspace preferred by ZHSEQR, as
      //   calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
      //   the worst case.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            MAXWRK = 1;
         } else {
            MAXWRK = N + N*ilaenv( 1, 'ZGEHRD', ' ', N, 1, N, 0 );
            MINWRK = 2*N;
            if ( WANTVL ) {
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ilaenv( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) );
               ztrevc3('L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               zhseqr('S', 'V', N, 1, N, A, LDA, W, VL, LDVL, WORK, -1, INFO );
            } else if ( WANTVR ) {
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ilaenv( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) );
               ztrevc3('R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR );
               LWORK_TREVC = INT( WORK(1) );
               MAXWRK = max( MAXWRK, N + LWORK_TREVC );
               zhseqr('S', 'V', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO );
            } else {
               zhseqr('E', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO );
            }
            HSWORK = INT( WORK(1) );
            MAXWRK = max( MAXWRK, HSWORK, MINWRK );
         }
         WORK[1] = MAXWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGEEV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = dlamch( 'P' );
      SMLNUM = dlamch( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM );
      SCALEA = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         SCALEA = true;
         CSCALE = SMLNUM;
      } else if ( ANRM > BIGNUM ) {
         SCALEA = true;
         CSCALE = BIGNUM;
      }
      if (SCALEA) zlascl( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR );

      // Balance the matrix
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1;
      zgebal('B', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR );

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1;
      IWRK = ITAU + N;
      zgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L';
         zlacpy('L', N, N, A, LDA, VL, LDVL );

         // Generate unitary matrix in VL
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         zunghr(N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VL
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         zhseqr('S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO );

         if ( WANTVR ) {

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B';
            zlacpy('F', N, N, VL, LDVL, VR, LDVR );
         }

      } else if ( WANTVR ) {

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R';
         zlacpy('L', N, N, A, LDA, VR, LDVR );

         // Generate unitary matrix in VR
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         zunghr(N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

         // Perform QR iteration, accumulating Schur vectors in VR
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         zhseqr('S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );

      } else {

         // Compute eigenvalues only
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU;
         zhseqr('E', 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO );
      }

      // If INFO != 0 from ZHSEQR, then quit

      if (INFO != 0) GO TO 50;

      if ( WANTVL || WANTVR ) {

         // Compute left and/or right eigenvectors
         // (CWorkspace: need 2*N, prefer N + 2*N*NB)
         // (RWorkspace: need 2*N)

         IRWORK = IBAL + N;
         ztrevc3(SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, RWORK( IRWORK ), N, IERR );
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors
         // (CWorkspace: none)
         // (RWorkspace: need N)

         zgebak('B', 'L', N, ILO, IHI, RWORK( IBAL ), N, VL, LDVL, IERR );

         // Normalize left eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 20
            SCL = ONE / DZNRM2( N, VL( 1, I ), 1 );
            zdscal(N, SCL, VL( 1, I ), 1 );
            for (K = 1; K <= N; K++) { // 10
               RWORK[IRWORK+K-1] = (VL( K, I )).toDouble()**2 + AIMAG( VL( K, I ) )**2;
            } // 10
            K = idamax( N, RWORK( IRWORK ), 1 );
            TMP = CONJG( VL( K, I ) ) / sqrt( RWORK( IRWORK+K-1 ) );
            zscal(N, TMP, VL( 1, I ), 1 );
            VL[K][I] = DCMPLX( (VL( K, I )).toDouble(), ZERO );
         } // 20
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors
         // (CWorkspace: none)
         // (RWorkspace: need N)

         zgebak('B', 'R', N, ILO, IHI, RWORK( IBAL ), N, VR, LDVR, IERR );

         // Normalize right eigenvectors and make largest component real

         for (I = 1; I <= N; I++) { // 40
            SCL = ONE / DZNRM2( N, VR( 1, I ), 1 );
            zdscal(N, SCL, VR( 1, I ), 1 );
            for (K = 1; K <= N; K++) { // 30
               RWORK[IRWORK+K-1] = (VR( K, I )).toDouble()**2 + AIMAG( VR( K, I ) )**2;
            } // 30
            K = idamax( N, RWORK( IRWORK ), 1 );
            TMP = CONJG( VR( K, I ) ) / sqrt( RWORK( IRWORK+K-1 ) );
            zscal(N, TMP, VR( 1, I ), 1 );
            VR[K][I] = DCMPLX( (VR( K, I )).toDouble(), ZERO );
         } // 40
      }

      // Undo scaling if necessary

      } // 50
      if ( SCALEA ) {
         zlascl('G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), max( N-INFO, 1 ), IERR );
         if ( INFO > 0 ) {
            zlascl('G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR );
         }
      }

      WORK[1] = MAXWRK;
      }
