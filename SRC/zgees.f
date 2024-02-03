      SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
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
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
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
        // HSWORK refers to the workspace preferred by ZHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N

            CALL ZHSEQR( 'S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL )
            HSWORK = INT( WORK( 1 ) )

            if ( .NOT.WANTVS ) {
               MAXWRK = MAX( MAXWRK, HSWORK )
            } else {
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, HSWORK )
            }
         }
         WORK( 1 ) = MAXWRK

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZGEES ', -INFO )
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

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      } else if ( ANRM.GT.BIGNUM ) {
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      }
      IF( SCALEA ) CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

      // Permute the matrix to make it more nearly triangular
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1
      CALL ZGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1
      IWRK = N + ITAU
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      if ( WANTVS ) {

         // Copy Householder vectors to VS

         CALL ZLACPY( 'L', N, N, A, LDA, VS, LDVS )

         // Generate unitary matrix in VS
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         CALL ZUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
      }

      SDIM = 0

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (CWorkspace: need 1, prefer HSWORK (see comments) )
      // (RWorkspace: none)

      IWRK = ITAU
      CALL ZHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL.GT.0 ) INFO = IEVAL

      // Sort eigenvalues if desired

      if ( WANTST .AND. INFO.EQ.0 ) {
         IF( SCALEA ) CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
         DO 10 I = 1, N
            BWORK( I ) = SELECT( W( I ) )
   10    CONTINUE

         // Reorder eigenvalues and transform Schur vectors
         // (CWorkspace: none)
         // (RWorkspace: none)

         CALL ZTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND )
      }

      if ( WANTVS ) {

         // Undo balancing
         // (CWorkspace: none)
         // (RWorkspace: need N)

         CALL ZGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS, IERR )
      }

      if ( SCALEA ) {

         // Undo scaling for the Schur form of A

         CALL ZLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL ZCOPY( N, A, LDA+1, W, 1 )
      }

      WORK( 1 ) = MAXWRK
      RETURN

      // End of ZGEES

      }
