      void cheev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * );
      COMPLEX            A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWKOPT, NB;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      //- REAL               CLANHE, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL ILAENV, LSAME, CLANHE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHETRD, CLASCL, CSTEQR, CUNGTR, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      LOWER = LSAME( UPLO, 'L' );
      LQUERY = ( LWORK == -1 );

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( LOWER || LSAME( UPLO, 'U' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }

      if ( INFO == 0 ) {
         NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 );
         LWKOPT = max( 1, ( NB+1 )*N );
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if( LWORK < max( 1, 2*N-1 ) && !LQUERY ) INFO = -8;
      }

      if ( INFO != 0 ) {
         xerbla('CHEEV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         return;
      }

      if ( N == 1 ) {
         W[1] = REAL( A( 1, 1 ) );
         WORK[1] = 1;
         if (WANTZ) A( 1, 1 ) = CONE;
         return;
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = CLANHE( 'M', UPLO, N, A, LDA, RWORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if (ISCALE == 1) clascl( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call CHETRD to reduce Hermitian matrix to tridiagonal form.

      INDE = 1;
      INDTAU = 1;
      INDWRK = INDTAU + N;
      LLWORK = LWORK - INDWRK + 1;
      chetrd(UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // CUNGTR to generate the unitary matrix, then call CSTEQR.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cungtr(UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
         INDWRK = INDE + N;
         csteqr(JOBZ, N, W, RWORK( INDE ), A, LDA, RWORK( INDWRK ), INFO );
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

      // Set WORK(1) to optimal complex workspace size.

      WORK[1] = SROUNDUP_LWORK(LWKOPT);

      return;
      }
