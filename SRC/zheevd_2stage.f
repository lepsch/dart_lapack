      SUBROUTINE ZHEEVD_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, INDWRK, ISCALE, LIWMIN, LLRWK, LLWORK, LLWRK2, LRWMIN, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS

       DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV2STAGE
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, DLAMCH, ZLANHE, ILAENV2STAGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZLACPY, ZLASCL, ZSTEDC, ZUNMTR, ZHETRD_2STAGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            LRWMIN = 1
            LIWMIN = 1
         ELSE
            KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
            IF( WANTZ ) THEN
               LWMIN = 2*N + N*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            ELSE
               LWMIN = N + 1 + LHTRD + LWTRD
               LRWMIN = N
               LIWMIN = 1
            END IF
         END IF
         WORK( 1 )  = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEVD_2STAGE', -INFO )
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
         W( 1 ) = DBLE( A( 1, 1 ) )
         IF( WANTZ ) A( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.
*
      INDE    = 1
      INDRWK  = INDE + N
      LLRWK   = LRWORK - INDRWK + 1
      INDTAU  = 1
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1
      INDWK2  = INDWRK + N*N
      LLWRK2  = LWORK - INDWK2 + 1
*
      CALL ZHETRD_2STAGE( JOBZ, UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
*     tridiagonal matrix, then call ZUNMTR to multiply it to the
*     Householder transformations represented as Householder vectors in
*     A.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZSTEDC( 'I', N, W, RWORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO )
         CALL ZUNMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL ZLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
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
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 )  = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of ZHEEVD_2STAGE
*
      END
