      void chbev_2stage(final int JOBZ, final int UPLO, final int N, final int KD, final Matrix<double> AB, final int LDAB, final int W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Box<int> INFO,) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N, LWORK;
      double               RWORK( * ), W( * );
      Complex            AB( LDAB, * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LOWER, WANTZ, LQUERY;
      int                IINFO, IMAX, INDE, INDWRK, INDRWK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS;
      double               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- REAL               SLAMCH, CLANHB, SROUNDUP_LWORK;
      // EXTERNAL lsame, SLAMCH, CLANHB, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTERF, XERBLA, CLASCL, CSTEQR, CHETRD_2STAGE, CHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, SQRT

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
      } else if ( KD < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD+1 ) {
         INFO = -6;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9;
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LWMIN = 1;
            WORK[1] = SROUNDUP_LWORK(LWMIN);
         } else {
            IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 );
            LWMIN = LHTRD + LWTRD;
            WORK[1] = SROUNDUP_LWORK(LWMIN);
         }

         if (LWORK < LWMIN && !LQUERY) INFO = -11;
      }

      if ( INFO != 0 ) {
         xerbla('CHBEV_2STAGE ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( LOWER ) {
            W[1] = double( AB( 1, 1 ) );
         } else {
            W[1] = double( AB( KD+1, 1 ) );
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
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

      // Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDE    = 1;
      INDHOUS = 1;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;

      chetrd_hb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

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

      // Set WORK(1) to optimal workspace size.

      WORK[1] = SROUNDUP_LWORK(LWMIN);

      }
