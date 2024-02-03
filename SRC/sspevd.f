      SUBROUTINE SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWRK, ISCALE, LIWMIN, LLWORK, LWMIN;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSP, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANSP, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SOPMTR, SSCAL, SSPTRD, SSTEDC, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -7
      }

      if ( INFO.EQ.0 ) {
         if ( N.LE.1 ) {
            LIWMIN = 1
            LWMIN = 1
         } else {
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N
               LWMIN = 1 + 6*N + N**2
            } else {
               LIWMIN = 1
               LWMIN = 2*N
            }
         }
         IWORK( 1 ) = LIWMIN
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -9
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSPEVD', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         W( 1 ) = AP( 1 )
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = SLANSP( 'M', UPLO, N, AP, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         CALL SSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      }

      // Call SSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDE = 1
      INDTAU = INDE + N
      CALL SSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call SOPMTR to multiply it by the
      // Householder transformations represented in AP.

      if ( .NOT.WANTZ ) {
         CALL SSTERF( N, W, WORK( INDE ), INFO )
      } else {
         INDWRK = INDTAU + N
         LLWORK = LWORK - INDWRK + 1
         CALL SSTEDC( 'I', N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWORK, IWORK, LIWORK, INFO )          CALL SOPMTR( 'L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO )
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      IF( ISCALE.EQ.1 ) CALL SSCAL( N, ONE / SIGMA, W, 1 )

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of SSPEVD

      }
