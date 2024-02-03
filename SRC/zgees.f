      void zgees(JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      double             RWORK( * );
      Complex         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * );
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTST, WANTVS;
      int                HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, ITAU, IWRK, MAXWRK, MINWRK;
      double             ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZTRSEN, ZUNGHR
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      WANTVS = LSAME( JOBVS, 'V' );
      WANTST = LSAME( SORT, 'S' );
      if ( ( !WANTVS ) && ( !LSAME( JOBVS, 'N' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTST ) && ( !LSAME( SORT, 'N' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDVS < 1 || ( WANTVS && LDVS < N ) ) {
         INFO = -10;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to real
        // workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by ZHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            MAXWRK = 1;
         } else {
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 );
            MINWRK = 2*N;

            zhseqr('S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL );
            HSWORK = INT( WORK( 1 ) );

            if ( !WANTVS ) {
               MAXWRK = max( MAXWRK, HSWORK );
            } else {
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) );
               MAXWRK = max( MAXWRK, HSWORK );
            }
         }
         WORK( 1 ) = MAXWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGEES ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         SDIM = 0;
         return;
      }

      // Get machine constants

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' );
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

      // Permute the matrix to make it more nearly triangular
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1;
      zgebal('P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR );

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1;
      IWRK = N + ITAU;
      zgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVS ) {

         // Copy Householder vectors to VS

         zlacpy('L', N, N, A, LDA, VS, LDVS );

         // Generate unitary matrix in VS
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         zunghr(N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );
      }

      SDIM = 0;

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (CWorkspace: need 1, prefer HSWORK (see comments) )
      // (RWorkspace: none)

      IWRK = ITAU;
      CALL ZHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL > 0 ) INFO = IEVAL;

      // Sort eigenvalues if desired

      if ( WANTST && INFO == 0 ) {
         if (SCALEA) zlascl( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR );
         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELECT( W( I ) );
         } // 10

         // Reorder eigenvalues and transform Schur vectors
         // (CWorkspace: none)
         // (RWorkspace: none)

         ztrsen('N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND );
      }

      if ( WANTVS ) {

         // Undo balancing
         // (CWorkspace: none)
         // (RWorkspace: need N)

         zgebak('P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS, IERR );
      }

      if ( SCALEA ) {

         // Undo scaling for the Schur form of A

         zlascl('U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR );
         zcopy(N, A, LDA+1, W, 1 );
      }

      WORK( 1 ) = MAXWRK;
      return;
      }
