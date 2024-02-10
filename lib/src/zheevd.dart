      void zheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double             RWORK( * ), W( * );
      Complex         A( LDA, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, INDWRK, ISCALE, LIOPT, LIWMIN, LLRWK, LLWORK, LLWRK2, LOPT, LROPT, LRWMIN, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANHE;
      // EXTERNAL lsame, ILAENV, DLAMCH, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZHETRD, ZLACPY, ZLASCL, ZSTEDC, ZUNMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( LOWER || lsame( UPLO, 'U' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LWMIN = 1;
            LRWMIN = 1;
            LIWMIN = 1;
            LOPT = LWMIN;
            LROPT = LRWMIN;
            LIOPT = LIWMIN;
         } else {
            if ( WANTZ ) {
               LWMIN = 2*N + N*N;
               LRWMIN = 1 + 5*N + 2*N**2;
               LIWMIN = 3 + 5*N;
            } else {
               LWMIN = N + 1;
               LRWMIN = N;
               LIWMIN = 1;
            }
            LOPT = max( LWMIN, N + N*ilaenv( 1, 'ZHETRD', UPLO, N, -1, -1, -1 ) );
            LROPT = LRWMIN;
            LIOPT = LIWMIN;
         }
         WORK[1] = LOPT;
         RWORK[1] = LROPT;
         IWORK[1] = LIOPT;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -10;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHEEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = (A( 1, 1 )).toDouble();
         if (WANTZ) A( 1, 1 ) = CONE;
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

      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if (ISCALE == 1) zlascl( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

      INDE = 1;
      INDTAU = 1;
      INDWRK = INDTAU + N;
      INDRWK = INDE + N;
      INDWK2 = INDWRK + N*N;
      LLWORK = LWORK - INDWRK + 1;
      LLWRK2 = LWORK - INDWK2 + 1;
      LLRWK = LRWORK - INDRWK + 1;
      zhetrd(UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call ZUNMTR to multiply it to the
      // Householder transformations represented as Householder vectors in
      // A.

      if ( !WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         zunmtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         zlacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = N;
         } else {
            IMAX = INFO - 1;
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      WORK[1] = LOPT;
      RWORK[1] = LROPT;
      IWORK[1] = LIOPT;

      }
