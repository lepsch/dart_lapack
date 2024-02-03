      SUBROUTINE DSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      int                IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLWRK2, LWMIN       DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DLAMCH, DLANSB
      EXTERNAL           LSAME, DLAMCH, DLANSB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASCL, DSBTRD, DSCAL, DSTEDC, DSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE
         IF( WANTZ ) THEN
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N
         END IF
      END IF
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = AB( 1, 1 )
         IF( WANTZ ) Z( 1, 1 ) = ONE
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
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = DLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            CALL DLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL DLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
      END IF
*
*     Call DSBTRD to reduce symmetric band matrix to tridiagonal form.
*
      INDE = 1
      INDWRK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1
      CALL DSBTRD( JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DSTEDC( 'I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )          CALL DGEMM( 'N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N )
         CALL DLACPY( 'A', N, N, WORK( INDWK2 ), N, Z, LDZ )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) CALL DSCAL( N, ONE / SIGMA, W, 1 )
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DSBEVD
*
      END
