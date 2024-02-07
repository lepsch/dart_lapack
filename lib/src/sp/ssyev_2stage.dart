      void ssyev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, LDA, LWORK, N;
      double               A( LDA, * ), W( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      double               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- REAL               SLAMCH, SLANSY, SROUNDUP_LWORK;
      // EXTERNAL lsame, SLAMCH, SLANSY, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASCL, SORGTR, SSCAL, SSTEQR, SSTERF, XERBLA, SSYTRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 );

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
         KD    = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1, -1 );
         IB    = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1, -1 );
         LHTRD = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 );
         LWTRD = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 );
         LWMIN = 2*N + LHTRD + LWTRD;
         WORK[1] = LWMIN;

         if (LWORK < LWMIN && !LQUERY) INFO = -8;
      }

      if ( INFO != 0 ) {
         xerbla('SSYEV_2STAGE ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         return;
      }

      if ( N == 1 ) {
         W[1] = A( 1, 1 );
         WORK[1] = 2;
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

      ssytrd_2stage(JOBZ, UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // SORGTR to generate the orthogonal matrix, then call SSTEQR.

      if ( !WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         // Not available in this release, and argument checking should not
         // let it getting here
         return;
         sorgtr(UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
         ssteqr(JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = N;
         } else {
            IMAX = INFO - 1;
         }
         sscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal workspace size.

      WORK[1] = SROUNDUP_LWORK(LWMIN);

      return;
      }
