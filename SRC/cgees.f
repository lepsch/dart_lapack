      SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTST, WANTVS;
      int                HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, ITAU, IWRK, MAXWRK, MINWRK;
      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY, CLASCL, CTRSEN, CUNGHR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      if ( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) {
         INFO = -1
      } else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) {
         INFO = -10
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to real
        // workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by CHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = N + N*ILAENV( 1, 'CGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N

            chseqr('S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL );
            HSWORK = INT( WORK( 1 ) )

            if ( .NOT.WANTVS ) {
               MAXWRK = MAX( MAXWRK, HSWORK )
            } else {
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'CUNGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, HSWORK )
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('CGEES ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         SDIM = 0
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      } else if ( ANRM.GT.BIGNUM ) {
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      }
      if (SCALEA) CALL CLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR );

      // Permute the matrix to make it more nearly triangular
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1
      cgebal('P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR );

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1
      IWRK = N + ITAU
      cgehrd(N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );

      if ( WANTVS ) {

         // Copy Householder vectors to VS

         clacpy('L', N, N, A, LDA, VS, LDVS );

         // Generate unitary matrix in VS
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         cunghr(N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR );
      }

      SDIM = 0

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (CWorkspace: need 1, prefer HSWORK (see comments) )
      // (RWorkspace: none)

      IWRK = ITAU
      CALL CHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL.GT.0 ) INFO = IEVAL

      // Sort eigenvalues if desired

      if ( WANTST .AND. INFO.EQ.0 ) {
         if (SCALEA) CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR );
         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELECT( W( I ) )
         } // 10

         // Reorder eigenvalues and transform Schur vectors
         // (CWorkspace: none)
         // (RWorkspace: none)

         ctrsen('N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND );
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
      }

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of CGEES

      }
