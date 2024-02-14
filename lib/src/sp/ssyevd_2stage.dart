      void ssyevd_2stage(final int JOBZ, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final int W, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final int LIWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LWORK, N;
      int                IWORK( * );
      double               A( LDA, * ), W( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;

      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE, LIWMIN, LLWORK, LLWRK2, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      double               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- REAL               SLAMCH, SLANSY;
      // EXTERNAL lsame, SLAMCH, SLANSY, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLASCL, SORMTR, SSCAL, SSTEDC, SSTERF, SSYTRD_2STAGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

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
            LIWMIN = 1;
            LWMIN = 1;
         } else {
            KD    = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 );
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N;
               LWMIN = 1 + 6*N + 2*N**2;
            } else {
               LIWMIN = 1;
               LWMIN = 2*N + 1 + LHTRD + LWTRD;
            }
         }
         WORK[1] = LWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -10;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSYEVD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = A( 1, 1 );
         if (WANTZ) A( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' );
      EPS    = SLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN   = sqrt( SMLNUM );
      RMAX   = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if (ISCALE == 1) slascl( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.

      INDE    = 1;
      INDTAU  = INDE + N;
      INDHOUS = INDTAU + N;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;
      INDWK2  = INDWRK + N*N;
      LLWRK2  = LWORK - INDWK2 + 1;

      ssytrd_2stage(JOBZ, UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call SORMTR to multiply it by the
      // Householder transformations stored in A.

      if ( !WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         // Not available in this release, and argument checking should not
         // let it getting here
         return;
         sstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         sormtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         slacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) sscal( N, ONE / SIGMA, W, 1 );

      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;

      }
