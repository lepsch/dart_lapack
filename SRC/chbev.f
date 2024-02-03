      void chbev(JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, RWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * );
      COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, ISCALE;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHB, SLAMCH;
      // EXTERNAL LSAME, CLANHB, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHBTRD, CLASCL, CSTEQR, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      LOWER = LSAME( UPLO, 'L' );

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( LOWER || LSAME( UPLO, 'U' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KD < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD+1 ) {
         INFO = -6;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9;
      }

      if ( INFO != 0 ) {
         xerbla('CHBEV ', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( LOWER ) {
            W( 1 ) = REAL( AB( 1, 1 ) );
         } else {
            W( 1 ) = REAL( AB( KD+1, 1 ) );
         }
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

      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK );
      ISCALE = 0;
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            clascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            clascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.

      INDE = 1;
      chbtrd(JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
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

      return;

      // End of CHBEV

      }
