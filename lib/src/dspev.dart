      void dspev(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANSP;
      // EXTERNAL lsame, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DOPGTR, DSCAL, DSPTRD, DSTEQR, DSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );

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

      if ( INFO != 0 ) {
         xerbla('DSPEV ', -INFO );
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

      SAFMIN = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = DLANSP( 'M', UPLO, N, AP, WORK );
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
      // DOPGTR to generate the orthogonal matrix, then call DSTEQR.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         INDWRK = INDTAU + N;
         dopgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
         dsteqr(JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDTAU ), INFO );
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

      return;
      }