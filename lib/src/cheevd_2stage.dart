      void cheevd_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) {

      // IMPLICIT NONE

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
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
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, INDWRK, ISCALE, LIWMIN, LLRWK, LLWORK, LLWRK2, LRWMIN, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;

       REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV2STAGE;
      //- REAL               SLAMCH, CLANHE;
      // EXTERNAL LSAME, SLAMCH, CLANHE, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTERF, XERBLA, CLACPY, CLASCL, CSTEDC, CUNMTR, CHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      LOWER = LSAME( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( !( LSAME( JOBZ, 'N' ) ) ) {
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
            LWMIN = 1;
            LRWMIN = 1;
            LIWMIN = 1;
         } else {
            KD    = ILAENV2STAGE( 1, 'CHETRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'CHETRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 );
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
         WORK( 1 )  = LWMIN;
         RWORK( 1 ) = LRWMIN;
         IWORK( 1 ) = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -10;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEEVD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W( 1 ) = REAL( A( 1, 1 ) );
         if (WANTZ) A( 1, 1 ) = CONE;
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

      // Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDE    = 1;
      INDRWK  = INDE + N;
      LLRWK   = LRWORK - INDRWK + 1;
      INDTAU  = 1;
      INDHOUS = INDTAU + N;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;
      INDWK2  = INDWRK + N*N;
      LLWRK2  = LWORK - INDWK2 + 1;

      chetrd_2stage(JOBZ, UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // CSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call CUNMTR to multiply it to the
      // Householder transformations represented as Householder vectors in
      // A.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cstedc('I', N, W, RWORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         cunmtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         clacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
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

      WORK( 1 )  = LWMIN;
      RWORK( 1 ) = LRWMIN;
      IWORK( 1 ) = LIWMIN;

      return;
      }
