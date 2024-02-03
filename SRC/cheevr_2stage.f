      SUBROUTINE CHEEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N;
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDIBL, INDIFL, INDISP, INDIWO, INDRD, INDRDD, INDRE, INDREE, INDRWK, INDTAU, INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, LLWORK, LLRWORK, LLWRKN, LRWMIN, LWMIN, NSPLIT, LHTRD, LWTRD, KD, IB, INDHOUS;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV, ILAENV2STAGE;
      REAL               SLAMCH, CLANSY, SROUNDUP_LWORK
      EXTERNAL           LSAME, SLAMCH, CLANSY, ILAENV, ILAENV2STAGE, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SSCAL, SSTEBZ, SSTERF, XERBLA, CSSCAL, CHETRD_2STAGE, CSTEMR, CSTEIN, CSWAP, CUNMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IEEEOK = ILAENV( 10, 'CHEEVR', 'N', 1, 2, 3, 4 )
*
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LRWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
*
      KD     = ILAENV2STAGE( 1, 'CHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'CHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
      LHTRD  = ILAENV2STAGE( 3, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
      LWTRD  = ILAENV2STAGE( 4, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
*
      IF( N.LE.1 ) THEN
         LWMIN  = 1
         LRWMIN = 1
         LIWMIN = 1
      ELSE
         LWMIN  = N + LHTRD + LWTRD
         LRWMIN = 24*N
         LIWMIN = 10*N
      END IF
*
      INFO = 0
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -8
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -15
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
         RWORK( 1 ) = SROUNDUP_LWORK( LRWMIN )
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -20
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -22
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHEEVR_2STAGE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         WORK( 1 ) = 1
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = REAL( A( 1, 1 ) )
         ELSE
            IF( VL.LT.REAL( A( 1, 1 ) ) .AND. VU.GE.REAL( A( 1, 1 ) ) ) THEN
               M = 1
               W( 1 ) = REAL( A( 1, 1 ) )
            END IF
         END IF
         IF( WANTZ ) THEN
            Z( 1, 1 ) = ONE
            ISUPPZ( 1 ) = 1
            ISUPPZ( 2 ) = 1
         END IF
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF (VALEIG) THEN
         VLL = VL
         VUU = VU
      END IF
      ANRM = CLANSY( 'M', UPLO, N, A, LDA, RWORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL CSSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL CSSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF

*     Initialize indices into workspaces.  Note: The IWORK indices are
*     used only if SSTERF or CSTEMR fail.

*     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
*     elementary reflectors used in CHETRD.
      INDTAU = 1
*     INDWK is the starting offset of the remaining complex workspace,
*     and LLWORK is the remaining complex workspace size.
      INDHOUS = INDTAU + N
      INDWK   = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWK + 1

*     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
*     entries.
      INDRD = 1
*     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
*     tridiagonal matrix from CHETRD.
      INDRE = INDRD + N
*     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
*     -written by CSTEMR (the SSTERF path copies the diagonal to W).
      INDRDD = INDRE + N
*     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
*     -written while computing the eigenvalues in SSTERF and CSTEMR.
      INDREE = INDRDD + N
*     INDRWK is the starting offset of the left-over real workspace, and
*     LLRWORK is the remaining workspace size.
      INDRWK = INDREE + N
      LLRWORK = LRWORK - INDRWK + 1

*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
*     stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
*     stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
*     that corresponding to eigenvectors that fail to converge in
*     CSTEIN.  This information is discarded; if any fail, the driver
*     returns INFO > 0.
      INDIFL = INDISP + N
*     INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N

*
*     Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.
*
      CALL CHETRD_2STAGE( JOBZ, UPLO, N, A, LDA, RWORK( INDRD ), RWORK( INDRE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO )
*
*     If all eigenvalues are desired
*     then call SSTERF or CSTEMR and CUNMTR.
*
      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.N ) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF( ( ALLEIG.OR.TEST ) .AND. ( IEEEOK.EQ.1 ) ) THEN
         IF( .NOT.WANTZ ) THEN
            CALL SCOPY( N, RWORK( INDRD ), 1, W, 1 )
            CALL SCOPY( N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 )
            CALL SSTERF( N, W, RWORK( INDREE ), INFO )
         ELSE
            CALL SCOPY( N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 )
            CALL SCOPY( N, RWORK( INDRD ), 1, RWORK( INDRDD ), 1 )
*
            IF ( ABSTOL .LE. TWO*N*EPS ) THEN
               TRYRAC = .TRUE.
            ELSE
               TRYRAC = .FALSE.
            END IF
            CALL CSTEMR( JOBZ, 'A', N, RWORK( INDRDD ), RWORK( INDREE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, RWORK( INDRWK ), LLRWORK, IWORK, LIWORK, INFO )
*
*           Apply unitary matrix used in reduction to tridiagonal
*           form to eigenvectors returned by CSTEMR.
*
            IF( WANTZ .AND. INFO.EQ.0 ) THEN
               INDWKN = INDWK
               LLWRKN = LWORK - INDWKN + 1
               CALL CUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
            END IF
         END IF
*
*
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
*     Also call SSTEBZ and CSTEIN if CSTEMR fails.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
       CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDRD ), RWORK( INDRE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL CSTEIN( N, RWORK( INDRD ), RWORK( INDRE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO )
*
*        Apply unitary matrix used in reduction to tridiagonal
*        form to eigenvectors returned by CSTEIN.
*
         INDWKN = INDWK
         LLWRKN = LWORK - INDWKN + 1
         CALL CUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
   30 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   40       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            END IF
   50    CONTINUE
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
      RWORK( 1 ) = SROUNDUP_LWORK( LRWMIN )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CHEEVR_2STAGE
*
      END
