      SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBVS, SENSE, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM
      DOUBLE PRECISION   RCONDE, RCONDV
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTSB, WANTSE, WANTSN, WANTST, WANTSV, WANTVS       int                HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, ITAU, IWRK, LWRK, MAXWRK, MINWRK
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL, XERBLA, ZCOPY, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZTRSEN, ZUNGHR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR. ( .NOT.WANTST .AND. .NOT.WANTSN ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
         INFO = -11
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of real workspace needed at that point in the
*       code, as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by ZHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.
*       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
*       depends on SDIM, which is computed by the routine ZTRSEN later
*       in the code.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            LWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N
*
            CALL ZHSEQR( 'S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS, WORK, -1, IEVAL )
            HSWORK = INT( WORK( 1 ) )
*
            IF( .NOT.WANTVS ) THEN
               MAXWRK = MAX( MAXWRK, HSWORK )
            ELSE
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, HSWORK )
            END IF
            LWRK = MAXWRK
            IF( .NOT.WANTSN ) LWRK = MAX( LWRK, ( N*N )/2 )
         END IF
         WORK( 1 ) = LWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -15
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEESX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
*
*     Permute the matrix to make it more nearly triangular
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      IBAL = 1
      CALL ZGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      ITAU = 1
      IWRK = N + ITAU
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
*
      IF( WANTVS ) THEN
*
*        Copy Householder vectors to VS
*
         CALL ZLACPY( 'L', N, N, A, LDA, VS, LDVS )
*
*        Generate unitary matrix in VS
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL ZUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
      END IF
*
      SDIM = 0
*
*     Perform QR iteration, accumulating Schur vectors in VS if desired
*     (CWorkspace: need 1, prefer HSWORK (see comments) )
*     (RWorkspace: none)
*
      IWRK = ITAU
      CALL ZHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL.GT.0 ) INFO = IEVAL
*
*     Sort eigenvalues if desired
*
      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA ) CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
         DO 10 I = 1, N
            BWORK( I ) = SELECT( W( I ) )
   10    CONTINUE
*
*        Reorder eigenvalues, transform Schur vectors, and compute
*        reciprocal condition numbers
*        (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
*                     otherwise, need none )
*        (RWorkspace: none)
*
         CALL ZTRSEN( SENSE, JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM, RCONDE, RCONDV, WORK( IWRK ), LWORK-IWRK+1, ICOND )
         IF( .NOT.WANTSN ) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
         IF( ICOND.EQ.-14 ) THEN
*
*           Not enough complex workspace
*
            INFO = -15
         END IF
      END IF
*
      IF( WANTVS ) THEN
*
*        Undo balancing
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL ZGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS, IERR )
      END IF
*
      IF( SCALEA ) THEN
*
*        Undo scaling for the Schur form of A
*
         CALL ZLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL ZCOPY( N, A, LDA+1, W, 1 )
         IF( ( WANTSV .OR. WANTSB ) .AND. INFO.EQ.0 ) THEN
            DUM( 1 ) = RCONDV
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
            RCONDV = DUM( 1 )
         END IF
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of ZGEESX
*
      END
