      void zgeesx(JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SENSE, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      double             RCONDE, RCONDV;
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
      bool               LQUERY, SCALEA, WANTSB, WANTSE, WANTSN, WANTST, WANTSV, WANTVS;
      int                HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, ITAU, IWRK, LWRK, MAXWRK, MINWRK;
      double             ANRM, BIGNUM, CSCALE, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, XERBLA, ZCOPY, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZTRSEN, ZUNGHR
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
      WANTVS = LSAME( JOBVS, 'V' );
      WANTST = LSAME( SORT, 'S' );
      WANTSN = LSAME( SENSE, 'N' );
      WANTSE = LSAME( SENSE, 'E' );
      WANTSV = LSAME( SENSE, 'V' );
      WANTSB = LSAME( SENSE, 'B' );
      LQUERY = ( LWORK == -1 );

      if ( ( !WANTVS ) && ( !LSAME( JOBVS, 'N' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTST ) && ( !LSAME( SORT, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( WANTSN || WANTSE || WANTSV || WANTSB ) || ( !WANTST && !WANTSN ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDVS < 1 || ( WANTVS && LDVS < N ) ) {
         INFO = -11;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of real workspace needed at that point in the
        // code, as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to real
        // workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by ZHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.
        // If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
        // depends on SDIM, which is computed by the routine ZTRSEN later
        // in the code.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            LWRK = 1;
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
            LWRK = MAXWRK;
            if ( !WANTSN) LWRK = max( LWRK, ( N*N )/2 );
         }
         WORK( 1 ) = LWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -15;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGEESX', -INFO );
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

         // Reorder eigenvalues, transform Schur vectors, and compute
         // reciprocal condition numbers
         // (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
                      // otherwise, need none )
         // (RWorkspace: none)

         ztrsen(SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, RCONDE, RCONDV, WORK( IWRK ), LWORK-IWRK+1, ICOND );
         if ( !WANTSN) MAXWRK = max( MAXWRK, 2*SDIM*( N-SDIM ) );
         if ( ICOND == -14 ) {

            // Not enough complex workspace

            INFO = -15;
         }
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
         if ( ( WANTSV || WANTSB ) && INFO == 0 ) {
            DUM( 1 ) = RCONDV;
            dlascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR );
            RCONDV = DUM( 1 );
         }
      }

      WORK( 1 ) = MAXWRK;
      return;
      }
