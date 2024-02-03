      SUBROUTINE CHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, N, LWORK
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ, LQUERY
      INTEGER            IINFO, IMAX, INDE, INDWRK, INDRWK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS       REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV2STAGE
      REAL               SLAMCH, CLANHB, SROUNDUP_LWORK
      EXTERNAL           LSAME, SLAMCH, CLANHB, ILAENV2STAGE, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSTERF, XERBLA, CLASCL, CSTEQR, CHETRD_2STAGE, CHETRD_HB2ST
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
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
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
         IF( N.LE.1 ) THEN
            LWMIN = 1
            WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         ELSE
            IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = LHTRD + LWTRD
            WORK( 1 )  = SROUNDUP_LWORK(LWMIN)
         ENDIF
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -11
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBEV_2STAGE ', -INFO )
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
         IF( LOWER ) THEN
            W( 1 ) = REAL( AB( 1, 1 ) )
         ELSE
            W( 1 ) = REAL( AB( KD+1, 1 ) )
         END IF
         IF( WANTZ ) Z( 1, 1 ) = ONE
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
      INDHOUS = 1
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1
*
      CALL CHETRD_HB2ST( "N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         INDRWK = INDE + N
         CALL CSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDRWK ), INFO )
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
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
*
      RETURN
*
*     End of CHBEV_2STAGE
*
      END
