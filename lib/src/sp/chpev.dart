      void chpev(final int JOBZ, final int UPLO, final int N, final int AP, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> _WORK_, final Array<double> RWORK_, final Box<int> INFO,) {
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, LDZ, N;
      double               RWORK( * ), W( * );
      Complex            AP( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK, ISCALE;
      double               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHP, SLAMCH;
      // EXTERNAL lsame, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPTRD, CSSCAL, CSTEQR, CUPGTR, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( lsame( UPLO, 'L' ) || lsame( UPLO, 'U' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         xerbla('CHPEV ', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = double( AP( 1 ) );
         RWORK[1] = 1;
         if (WANTZ) Z( 1, 1 ) = ONE;
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

      ANRM = CLANHP( 'M', UPLO, N, AP, RWORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         csscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form.

      INDE = 1;
      INDTAU = 1;
      chptrd(UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // CUPGTR to generate the orthogonal matrix, then call CSTEQR.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         INDWRK = INDTAU + N;
         cupgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
         INDRWK = INDE + N;
         csteqr(JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDRWK ), INFO );
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

      }
