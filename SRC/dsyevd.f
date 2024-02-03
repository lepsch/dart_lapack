
// >
// =====================================================================
      void dsyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), W( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..

      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE, LIOPT, LIWMIN, LLWORK, LLWRK2, LOPT, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DLASCL, DORMTR, DSCAL, DSTEDC, DSTERF, DSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      LOWER = LSAME( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

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
         if ( N <= 1 ) {
            LIWMIN = 1;
            LWMIN = 1;
            LOPT = LWMIN;
            LIOPT = LIWMIN;
         } else {
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N;
               LWMIN = 1 + 6*N + 2*N**2;
            } else {
               LIWMIN = 1;
               LWMIN = 2*N + 1;
            }
            LOPT = max( LWMIN, 2*N + N*ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ) );
            LIOPT = LIWMIN;
         }
         WORK( 1 ) = LOPT;
         IWORK( 1 ) = LIOPT;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -10;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSYEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W( 1 ) = A( 1, 1 );
         if (WANTZ) A( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if (ISCALE == 1) CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call DSYTRD to reduce symmetric matrix to tridiagonal form.

      INDE = 1;
      INDTAU = INDE + N;
      INDWRK = INDTAU + N;
      LLWORK = LWORK - INDWRK + 1;
      INDWK2 = INDWRK + N*N;
      LLWRK2 = LWORK - INDWK2 + 1;

      dsytrd(UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call DORMTR to multiply it by the
      // Householder transformations stored in A.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         dormtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         dlacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL DSCAL( N, ONE / SIGMA, W, 1 );

      WORK( 1 ) = LOPT;
      IWORK( 1 ) = LIOPT;

      return;

      // End of DSYEVD

      }
