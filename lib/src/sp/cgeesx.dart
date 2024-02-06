      void cgeesx(JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVS, SENSE, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      double               RCONDE, RCONDV;
      bool               BWORK( * );
      double               RWORK( * );
      Complex            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * );
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, SCALEA, WANTSB, WANTSE, WANTSN, WANTST, WANTSV, WANTVS;
      int                HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, ITAU, IWRK, LWRK, MAXWRK, MINWRK;
      double               ANRM, BIGNUM, CSCALE, EPS, SMLNUM;
      double               DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY, CLASCL, CTRSEN, CUNGHR, SLASCL, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               CLANGE, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input arguments

      INFO = 0;
      WANTVS = lsame( JOBVS, 'V' );
      WANTST = lsame( SORT, 'S' );
      WANTSN = lsame( SENSE, 'N' );
      WANTSE = lsame( SENSE, 'E' );
      WANTSV = lsame( SENSE, 'V' );
      WANTSB = lsame( SENSE, 'B' );
      LQUERY = ( LWORK == -1 );

      if ( ( !WANTVS ) && ( !lsame( JOBVS, 'N' ) ) ) {
         INFO = -1;
      } else if ( ( !WANTST ) && ( !lsame( SORT, 'N' ) ) ) {
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
        // HSWORK refers to the workspace preferred by CHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.
        // If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
        // depends on SDIM, which is computed by the routine CTRSEN later
        // in the code.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1;
            LWRK = 1;
         } else {
            MAXWRK = N + N*ilaenv( 1, 'CGEHRD', ' ', N, 1, N, 0 );
            MINWRK = 2*N;

            chseqr('S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL );
            HSWORK = INT( WORK( 1 ) );

            if ( !WANTVS ) {
               MAXWRK = max( MAXWRK, HSWORK );
            } else {
               MAXWRK = max( MAXWRK, N + ( N - 1 )*ilaenv( 1, 'CUNGHR', ' ', N, 1, N, -1 ) );
               MAXWRK = max( MAXWRK, HSWORK );
            }
            LWRK = MAXWRK;
            if ( !WANTSN) LWRK = max( LWRK, ( N*N )/2 );
         }
         WORK[1] = SROUNDUP_LWORK(LWRK);

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -15;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGEESX', -INFO );
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

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

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


      // Permute the matrix to make it more nearly triangular
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1;
      cgebal('P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR );

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1;
      IWRK = N + ITAU;
      cgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVS ) {

         // Copy Householder vectors to VS

         clacpy('L', N, N, A, LDA, VS, LDVS );

         // Generate unitary matrix in VS
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         cunghr(N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );
      }

      SDIM = 0;

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (CWorkspace: need 1, prefer HSWORK (see comments) )
      // (RWorkspace: none)

      IWRK = ITAU;
      CALL CHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL > 0 ) INFO = IEVAL;

      // Sort eigenvalues if desired

      if ( WANTST && INFO == 0 ) {
         if (SCALEA) clascl( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR );
         for (I = 1; I <= N; I++) { // 10
            BWORK[I] = SELECT( W( I ) );
         } // 10

         // Reorder eigenvalues, transform Schur vectors, and compute
         // reciprocal condition numbers
         // (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
                      // otherwise, need none )
         // (RWorkspace: none)

         ctrsen(SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, RCONDE, RCONDV, WORK( IWRK ), LWORK-IWRK+1, ICOND );
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

         cgebak('P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS, IERR );
      }

      if ( SCALEA ) {

         // Undo scaling for the Schur form of A

         clascl('U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR );
         ccopy(N, A, LDA+1, W, 1 );
         if ( ( WANTSV || WANTSB ) && INFO == 0 ) {
            DUM[1] = RCONDV;
            slascl('G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR );
            RCONDV = DUM( 1 );
         }
      }

      WORK[1] = SROUNDUP_LWORK(MAXWRK);
      return;
      }
