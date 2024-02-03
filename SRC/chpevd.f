      SUBROUTINE CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LRWORK, LWORK, N;
*     ..
*     .. Array Arguments ..
      int                IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            AP( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK, ISCALE, LIWMIN, LLRWK, LLWRK, LRWMIN, LWMIN       REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SLAMCH, SROUNDUP_LWORK
      EXTERNAL           LSAME, CLANHP, SLAMCH, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHPTRD, CSSCAL, CSTEDC, CUPMTR, SSCAL, SSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         ELSE
            IF( WANTZ ) THEN
               LWMIN = 2*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            ELSE
               LWMIN = N
               LRWMIN = N
               LIWMIN = 1
            END IF
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -9
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPEVD', -INFO )
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
         W( 1 ) = REAL( AP( 1 ) )
         IF( WANTZ ) Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = CLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL CSSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      END IF
*
*     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      INDRWK = INDE + N
      INDWRK = INDTAU + N
      LLWRK = LWORK - INDWRK + 1
      LLRWK = LRWORK - INDRWK + 1
      CALL CHPTRD( UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ), IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     CUPGTR to generate the orthogonal matrix, then call CSTEDC.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEDC( 'I', N, W, RWORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWRK, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO )
         CALL CUPMTR( 'L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO )
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
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of CHPEVD
*
      END
