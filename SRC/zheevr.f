      SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC
      CHARACTER          ORDER
      int                I, IEEEOK, IINFO, IMAX, INDIBL, INDIFL, INDISP, INDIWO, INDRD, INDRDD, INDRE, INDREE, INDRWK, INDTAU, INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, LLWORK, LLRWORK, LLWRKN, LRWMIN, LWKOPT, LWMIN, NB, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZHETRD, ZSTEMR, ZSTEIN, ZSWAP, ZUNMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IEEEOK = ILAENV( 10, 'ZHEEVR', 'N', 1, 2, 3, 4 )
*
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LRWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
*
      IF( N.LE.1 ) THEN
         LWMIN  = 1
         LRWMIN = 1
         LIWMIN = 1
      ELSE
         LWMIN  = 2*N
         LRWMIN = 24*N
         LIWMIN = 10*N
      END IF
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
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
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         NB = MAX( NB, ILAENV( 1, 'ZUNMTR', UPLO, N, -1, -1, -1 ) )
         LWKOPT = MAX( ( NB+1 )*N, LWMIN )
         WORK( 1 )  = LWKOPT
         RWORK( 1 ) = LRWMIN
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
         CALL XERBLA( 'ZHEEVR', -INFO )
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
            W( 1 ) = DBLE( A( 1, 1 ) )
         ELSE
            IF( VL.LT.DBLE( A( 1, 1 ) ) .AND. VU.GE.DBLE( A( 1, 1 ) ) ) THEN
               M = 1
               W( 1 ) = DBLE( A( 1, 1 ) )
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
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF (VALEIG) THEN
         VLL = VL
         VUU = VU
      END IF
      ANRM = ZLANSY( 'M', UPLO, N, A, LDA, RWORK )
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
               CALL ZDSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL ZDSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF

*     Initialize indices into workspaces.  Note: The IWORK indices are
*     used only if DSTERF or ZSTEMR fail.

*     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
*     elementary reflectors used in ZHETRD.
      INDTAU = 1
*     INDWK is the starting offset of the remaining complex workspace,
*     and LLWORK is the remaining complex workspace size.
      INDWK = INDTAU + N
      LLWORK = LWORK - INDWK + 1

*     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
*     entries.
      INDRD = 1
*     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
*     tridiagonal matrix from ZHETRD.
      INDRE = INDRD + N
*     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
*     -written by ZSTEMR (the DSTERF path copies the diagonal to W).
      INDRDD = INDRE + N
*     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
*     -written while computing the eigenvalues in DSTERF and ZSTEMR.
      INDREE = INDRDD + N
*     INDRWK is the starting offset of the left-over real workspace, and
*     LLRWORK is the remaining workspace size.
      INDRWK = INDREE + N
      LLRWORK = LRWORK - INDRWK + 1

*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
*     stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
*     stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
*     that corresponding to eigenvectors that fail to converge in
*     DSTEIN.  This information is discarded; if any fail, the driver
*     returns INFO > 0.
      INDIFL = INDISP + N
*     INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N

*
*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
*
      CALL ZHETRD( UPLO, N, A, LDA, RWORK( INDRD ), RWORK( INDRE ), WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO )
*
*     If all eigenvalues are desired
*     then call DSTERF or ZSTEMR and ZUNMTR.
*
      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.N ) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF( ( ALLEIG.OR.TEST ) .AND. ( IEEEOK.EQ.1 ) ) THEN
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N, RWORK( INDRD ), 1, W, 1 )
            CALL DCOPY( N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 )
            CALL DSTERF( N, W, RWORK( INDREE ), INFO )
         ELSE
            CALL DCOPY( N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 )
            CALL DCOPY( N, RWORK( INDRD ), 1, RWORK( INDRDD ), 1 )
*
            IF (ABSTOL .LE. TWO*N*EPS) THEN
               TRYRAC = .TRUE.
            ELSE
               TRYRAC = .FALSE.
            END IF
            CALL ZSTEMR( JOBZ, 'A', N, RWORK( INDRDD ), RWORK( INDREE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, RWORK( INDRWK ), LLRWORK, IWORK, LIWORK, INFO )
*
*           Apply unitary matrix used in reduction to tridiagonal
*           form to eigenvectors returned by ZSTEMR.
*
            IF( WANTZ .AND. INFO.EQ.0 ) THEN
               INDWKN = INDWK
               LLWRKN = LWORK - INDWKN + 1
               CALL ZUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
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
*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
*     Also call DSTEBZ and ZSTEIN if ZSTEMR fails.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
       CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDRD ), RWORK( INDRE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL ZSTEIN( N, RWORK( INDRD ), RWORK( INDRE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO )
*
*        Apply unitary matrix used in reduction to tridiagonal
*        form to eigenvectors returned by ZSTEIN.
*
         INDWKN = INDWK
         LLWRKN = LWORK - INDWKN + 1
         CALL ZUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
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
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
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
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            END IF
   50    CONTINUE
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 )  = LWKOPT
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of ZHEEVR
*
      END
