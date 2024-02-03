      SUBROUTINE SSYEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * )
      REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ, TRYRAC, TEST
      CHARACTER          ORDER
      int                I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE, INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU, INDWK, INDWKN, ISCALE, J, JJ, LIWMIN, LLWORK, LLWRKN, LWMIN, NSPLIT, LHTRD, LWTRD, KD, IB, INDHOUS
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV, ILAENV2STAGE
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      EXTERNAL           LSAME, SLAMCH, SLANSY, SROUNDUP_LWORK, ILAENV, ILAENV2STAGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SORMTR, SSCAL, SSTEBZ, SSTEMR, SSTEIN, SSTERF, SSWAP, SSYTRD_2STAGE, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IEEEOK = ILAENV( 10, 'SSYEVR', 'N', 1, 2, 3, 4 )
*
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
*
      KD     = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )
      LHTRD  = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
      LWTRD  = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
*
      IF( N.LE.1 ) THEN
         LWMIN  = 1
         LIWMIN = 1
      ELSE
         LWMIN  = MAX( 26*N, 5*N + LHTRD + LWTRD )
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
         ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -20
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
*         NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
*         NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) )
*         LWKOPT = MAX( ( NB+1 )*N, LWMIN )
         WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYEVR_2STAGE', -INFO )
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
            W( 1 ) = A( 1, 1 )
         ELSE
            IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
               M = 1
               W( 1 ) = A( 1, 1 )
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
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
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
               CALL SSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL SSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF

*     Initialize indices into workspaces.  Note: The IWORK indices are
*     used only if SSTERF or SSTEMR fail.

*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
*     elementary reflectors used in SSYTRD.
      INDTAU = 1
*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
      INDD = INDTAU + N
*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
*     tridiagonal matrix from SSYTRD.
      INDE = INDD + N
*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
*     -written by SSTEMR (the SSTERF path copies the diagonal to W).
      INDDD = INDE + N
*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
*     -written while computing the eigenvalues in SSTERF and SSTEMR.
      INDEE = INDDD + N
*     INDHOUS is the starting offset Householder storage of stage 2
      INDHOUS = INDEE + N
*     INDWK is the starting offset of the left-over workspace, and
*     LLWORK is the remaining workspace size.
      INDWK  = INDHOUS + LHTRD
      LLWORK = LWORK - INDWK + 1


*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
*     stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
*     stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
*     that corresponding to eigenvectors that fail to converge in
*     SSTEIN.  This information is discarded; if any fail, the driver
*     returns INFO > 0.
      INDIFL = INDISP + N
*     INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N

*
*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.
*
*
      CALL SSYTRD_2STAGE( JOBZ, UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO )
*
*     If all eigenvalues are desired
*     then call SSTERF or SSTEMR and SORMTR.
*
      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.N ) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF( ( ALLEIG.OR.TEST ) .AND. ( IEEEOK.EQ.1 ) ) THEN
         IF( .NOT.WANTZ ) THEN
            CALL SCOPY( N, WORK( INDD ), 1, W, 1 )
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SCOPY( N, WORK( INDD ), 1, WORK( INDDD ), 1 )
*
            IF (ABSTOL .LE. TWO*N*EPS) THEN
               TRYRAC = .TRUE.
            ELSE
               TRYRAC = .FALSE.
            END IF
            CALL SSTEMR( JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( INDWK ), LWORK, IWORK, LIWORK, INFO )
*
*
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by SSTEMR.
*
            IF( WANTZ .AND. INFO.EQ.0 ) THEN
               INDWKN = INDE
               LLWRKN = LWORK - INDWKN + 1
               CALL SORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
            END IF
         END IF
*
*
         IF( INFO.EQ.0 ) THEN
*           Everything worked.  Skip SSTEBZ/SSTEIN.  IWORK(:) are
*           undefined.
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.
*     Also call SSTEBZ and SSTEIN if SSTEMR fails.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
       CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ), IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL SSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO )
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by SSTEIN.
*
         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL SORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
*  Jump here if SSTEMR/SSTEIN succeeded.
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
*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
*     It may not be initialized (if SSTEMR/SSTEIN succeeded), and we do
*     not return this detailed information to the user.
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
               W( I ) = W( J )
               W( J ) = TMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            END IF
   50    CONTINUE
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of SSYEVR_2STAGE
*
      END
