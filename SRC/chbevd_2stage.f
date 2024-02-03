      SUBROUTINE CHBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDWK2, INDRWK, ISCALE, LLWORK, INDWK, LHTRD, LWTRD, IB, INDHOUS, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV2STAGE
      REAL               SLAMCH, CLANHB
      EXTERNAL           LSAME, SLAMCH, CLANHB, ILAENV2STAGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSTERF, XERBLA, CGEMM, CLACPY, CLASCL, CSTEDC, CHETRD_HB2ST
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LWMIN = 1
         LRWMIN = 1
         LIWMIN = 1
      ELSE
         IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
         IF( WANTZ ) THEN
            LWMIN = 2*N**2
            LRWMIN = 1 + 5*N + 2*N**2
            LIWMIN = 3 + 5*N
         ELSE
            LWMIN  = MAX( N, LHTRD + LWTRD )
            LRWMIN = N
            LIWMIN = 1
         END IF
      END IF
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
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
         WORK( 1 )  = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBEVD_2STAGE', -INFO )
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
         W( 1 ) = REAL( AB( 1, 1 ) )
         IF( WANTZ ) Z( 1, 1 ) = CONE
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
      RMAX   = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
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
            CALL CLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL CLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
      END IF
*
*     Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.
*
      INDE    = 1
      INDRWK  = INDE + N
      LLRWK   = LRWORK - INDRWK + 1
      INDHOUS = 1
      INDWK   = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWK + 1
      INDWK2  = INDWK + N*N
      LLWK2   = LWORK - INDWK2 + 1
*
      CALL CHETRD_HB2ST( "N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEDC( 'I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO )
         CALL CGEMM( 'N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N )
         CALL CLACPY( 'A', N, N, WORK( INDWK2 ), N, Z, LDZ )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 )  = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of CHBEVD_2STAGE
*
      END
