      void dsbevd_2stage(JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO ) {

      // IMPLICIT NONE

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, LLWRK2;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- double             DLAMCH, DLANSB;
      // EXTERNAL lsame, DLAMCH, DLANSB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASCL, DSCAL, DSTEDC, DSTERF, XERBLA, DSYTRD_SB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      LOWER = lsame( UPLO, 'L' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LIWMIN = 1;
         LWMIN = 1;
      } else {
         IB    = ILAENV2STAGE( 2, 'DSYTRD_SB2ST', JOBZ, N, KD, -1, -1 );
         LHTRD = ILAENV2STAGE( 3, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1 );
         LWTRD = ILAENV2STAGE( 4, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1 );
         if ( WANTZ ) {
            LIWMIN = 3 + 5*N;
            LWMIN = 1 + 5*N + 2*N**2;
         } else {
            LIWMIN = 1;
            LWMIN = max( 2*N, N+LHTRD+LWTRD );
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
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -13;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSBEVD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         W[1] = AB( 1, 1 );
         if (WANTZ) Z( 1, 1 ) = ONE;
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

      ANRM = DLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK );
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
            dlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            dlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call DSYTRD_SB2ST to reduce band symmetric matrix to tridiagonal form.

      INDE    = 1;
      INDHOUS = INDE + N;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;
      INDWK2  = INDWRK + N*N;
      LLWRK2  = LWORK - INDWK2 + 1;

      dsytrd_sb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         dgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N );
         dlacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) dscal( N, ONE / SIGMA, W, 1 );

      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;
      return;
      }
