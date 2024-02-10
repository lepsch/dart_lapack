      void zheevd_2stage(final int JOBZ, final int UPLO, final int N, final Matrix<double> A, final int LDA, final int W, final Array<double> WORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Array<int> IWORK, final int LIWORK, final Box<int> INFO ) {

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
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, INDWRK, ISCALE, LIWMIN, LLRWK, LLWORK, LLWRK2, LRWMIN, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;

       double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- double             DLAMCH, ZLANHE;
      // EXTERNAL lsame, DLAMCH, ZLANHE, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZLACPY, ZLASCL, ZSTEDC, ZUNMTR, ZHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( !( lsame( JOBZ, 'N' ) ) ) {
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
         } else {
            KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 );
            if ( WANTZ ) {
               LWMIN = 2*N + N*N;
               LRWMIN = 1 + 5*N + 2*N**2;
               LIWMIN = 3 + 5*N;
            } else {
               LWMIN = N + 1 + LHTRD + LWTRD;
               LRWMIN = N;
               LIWMIN = 1;
            }
         }
         WORK[1] = LWMIN;
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -10;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHEEVD_2STAGE', -INFO );
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
      EPS    = dlamch( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN   = sqrt( SMLNUM );
      RMAX   = sqrt( BIGNUM );

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

      // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDE    = 1;
      INDRWK  = INDE + N;
      LLRWK   = LRWORK - INDRWK + 1;
      INDTAU  = 1;
      INDHOUS = INDTAU + N;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;
      INDWK2  = INDWRK + N*N;
      LLWRK2  = LWORK - INDWK2 + 1;

      zhetrd_2stage(JOBZ, UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

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

      WORK[1] = LWMIN;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;

      }
