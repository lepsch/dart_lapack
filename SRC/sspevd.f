      SUBROUTINE SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWRK, ISCALE, LIWMIN, LLWORK, LWMIN;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSP, SROUNDUP_LWORK;
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

      WANTZ = LSAME( JOBZ, 'V' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( LSAME( UPLO, 'U' ) || LSAME( UPLO, 'L' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -7;
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LIWMIN = 1;
            LWMIN = 1;
         } else {
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N;
               LWMIN = 1 + 6*N + N**2;
            } else {
               LIWMIN = 1;
               LWMIN = 2*N;
            }
         }
         IWORK( 1 ) = LIWMIN;
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN);

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -9;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -11;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSPEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W( 1 ) = AP( 1 );
         if (WANTZ) Z( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = SQRT( SMLNUM );
      RMAX = SQRT( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = SLANSP( 'M', UPLO, N, AP, WORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         sscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call SSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDE = 1;
      INDTAU = INDE + N;
      ssptrd(UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call SOPMTR to multiply it by the
      // Householder transformations represented in AP.

      if ( !WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         INDWRK = INDTAU + N;
         LLWORK = LWORK - INDWRK + 1;
         sstedc('I', N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWORK, IWORK, LIWORK, INFO );
         sopmtr('L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL SSCAL( N, ONE / SIGMA, W, 1 );

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN);
      IWORK( 1 ) = LIWMIN;
      return;

      // End of SSPEVD

      }
