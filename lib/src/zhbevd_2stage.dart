      void zhbevd_2stage(JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double             RWORK( * ), W( * );
      Complex         AB( LDAB, * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDWK2, INDRWK, ISCALE, LLWORK, INDWK, LHTRD, LWTRD, IB, INDHOUS, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- double             DLAMCH, ZLANHB;
      // EXTERNAL lsame, DLAMCH, ZLANHB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZGEMM, ZLACPY, ZLASCL, ZSTEDC, ZHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, SQRT

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
         IB    = ILAENV2STAGE( 2, 'ZHETRD_HB2ST', JOBZ, N, KD, -1, -1 );
         LHTRD = ILAENV2STAGE( 3, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 );
         LWTRD = ILAENV2STAGE( 4, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 );
         if ( WANTZ ) {
            LWMIN = 2*N**2;
            LRWMIN = 1 + 5*N + 2*N**2;
            LIWMIN = 3 + 5*N;
         } else {
            LWMIN  = max( N, LHTRD + LWTRD );
            LRWMIN = N;
            LIWMIN = 1;
         }
      }
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
         WORK[1] = LWMIN;
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
         xerbla('ZHBEVD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = (AB( 1, 1 )).toDouble();
         if (WANTZ) Z( 1, 1 ) = CONE;
         return;
      }

      // Get machine constants.

      SAFMIN = dlamch( 'Safe minimum' );
      EPS    = dlamch( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN   = sqrt( SMLNUM );
      RMAX   = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ANRM = ZLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK );
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
            zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDE    = 1;
      INDRWK  = INDE + N;
      LLRWK   = LRWORK - INDRWK + 1;
      INDHOUS = 1;
      INDWK   = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWK + 1;
      INDWK2  = INDWK + N*N;
      LLWK2   = LWORK - INDWK2 + 1;

      zhetrd_hb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.

      if ( !WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         zgemm('N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N );
         zlacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
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

      WORK[1] = LWMIN;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;
      return;
      }
