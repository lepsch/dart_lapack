      SUBROUTINE SSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..

      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE, LIOPT, LIWMIN, LLWORK, LLWRK2, LOPT, LWMIN;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SLAMCH, SLANSY, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLASCL, SORMTR, SSCAL, SSTEDC, SSTERF, SSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER || LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LIWMIN = 1
            LWMIN = 1
            LOPT = LWMIN
            LIOPT = LIWMIN
         } else {
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N
               LWMIN = 1 + 6*N + 2*N**2
            } else {
               LIWMIN = 1
               LWMIN = 2*N + 1
            }
            LOPT = MAX( LWMIN, 2*N + N*ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 ) )
            LIOPT = LIWMIN
         }
         WORK( 1 ) = SROUNDUP_LWORK( LOPT )
         IWORK( 1 ) = LIOPT

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -8
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -10
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSYEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         W( 1 ) = A( 1, 1 )
         if (WANTZ) A( 1, 1 ) = ONE;
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

      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM > RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if (ISCALE == 1) CALL SLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call SSYTRD to reduce symmetric matrix to tridiagonal form.

      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1

      ssytrd(UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call SORMTR to multiply it by the
      // Householder transformations stored in A.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         sstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         sormtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         slacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL SSCAL( N, ONE / SIGMA, W, 1 );

      WORK( 1 ) = SROUNDUP_LWORK( LOPT )
      IWORK( 1 ) = LIOPT

      RETURN

      // End of SSYEVD

      }
