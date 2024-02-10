      void dspevd(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LWORK, N;
      int                IWORK( * );
      double             AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWRK, ISCALE, LIWMIN, LLWORK, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANSP;
      // EXTERNAL lsame, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DOPMTR, DSCAL, DSPTRD, DSTEDC, DSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( lsame( UPLO, 'U' ) || lsame( UPLO, 'L' ) ) ) {
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
         IWORK[1] = LIWMIN;
         WORK[1] = LWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -9;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -11;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSPEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = AP( 1 );
         if (WANTZ) Z( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = dlamch( 'Safe minimum' );
      EPS = dlamch( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = dlansp( 'M', UPLO, N, AP, WORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         dscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDE = 1;
      INDTAU = INDE + N;
      dsptrd(UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call DOPMTR to multiply it by the
      // Householder transformations represented in AP.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         INDWRK = INDTAU + N;
         LLWORK = LWORK - INDWRK + 1;
         dstedc('I', N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWORK, IWORK, LIWORK, INFO );
         dopmtr('L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) dscal( N, ONE / SIGMA, W, 1 );

      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;
      }
