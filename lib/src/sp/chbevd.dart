      void chbevd(JOBZ, UPLO, N, KD, final Matrix<double> AB, final int LDAB, W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Array<int> IWORK, final int LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double               RWORK( * ), W( * );
      Complex            AB( LDAB, * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
      double               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHB, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL lsame, CLANHB, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHBTRD, CLACPY, CLASCL, CSTEDC, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 || LRWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LWMIN = 1;
         LRWMIN = 1;
         LIWMIN = 1;
      } else {
         if ( WANTZ ) {
            LWMIN = 2*N**2;
            LRWMIN = 1 + 5*N + 2*N**2;
            LIWMIN = 3 + 5*N;
         } else {
            LWMIN = N;
            LRWMIN = N;
            LIWMIN = 1;
         }
      }
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
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
         WORK[1] = SROUNDUP_LWORK(LWMIN);
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -13;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -15;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHBEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = double( AB( 1, 1 ) );
         if (WANTZ) Z( 1, 1 ) = CONE;
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
      INDWRK = INDE + N;
      INDWK2 = 1 + N*N;
      LLWK2 = LWORK - INDWK2 + 1;
      LLRWK = LRWORK - INDWRK + 1;
      chbtrd(JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cstedc('I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDWRK ), LLRWK, IWORK, LIWORK, INFO );
         cgemm('N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N );
         clacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
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

      WORK[1] = SROUNDUP_LWORK(LWMIN);
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;
      }
