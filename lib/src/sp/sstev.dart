      void sstev(final int JOBZ, final int N, final int D, final int E, final Matrix<double> Z, final int LDZ, final Array<double> _WORK, final Box<int> INFO,) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ;
      int                INFO, LDZ, N;
      double               D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               WANTZ;
      int                IMAX, ISCALE;
      double               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TNRM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANST;
      // EXTERNAL lsame, SLAMCH, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTEQR, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -6;
      }

      if ( INFO != 0 ) {
         xerbla('SSTEV ', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
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

      ISCALE = 0;
      TNRM = SLANST( 'M', N, D, E );
      if ( TNRM > ZERO && TNRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / TNRM;
      } else if ( TNRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / TNRM;
      }
      if ( ISCALE == 1 ) {
         sscal(N, SIGMA, D, 1 );
         sscal(N-1, SIGMA, E( 1 ), 1 );
      }

      // For eigenvalues only, call SSTERF.  For eigenvalues and
      // eigenvectors, call SSTEQR.

      if ( !WANTZ ) {
         ssterf(N, D, E, INFO );
      } else {
         ssteqr('I', N, D, E, Z, LDZ, WORK, INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = N;
         } else {
            IMAX = INFO - 1;
         }
         sscal(IMAX, ONE / SIGMA, D, 1 );
      }

      }
