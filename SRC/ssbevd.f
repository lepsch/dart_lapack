      SUBROUTINE SSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLWRK2, LWMIN;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSB, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANSB, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASCL, SSBTRD, SSCAL, SSTEDC, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( N <= 1 ) {
         LIWMIN = 1
         LWMIN = 1
      } else {
         if ( WANTZ ) {
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N**2
         } else {
            LIWMIN = 1
            LWMIN = 2*N
         }
      }
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER || LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( KD < 0 ) {
         INFO = -4
      } else if ( LDAB < KD+1 ) {
         INFO = -6
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSBEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         W( 1 ) = AB( 1, 1 )
         if (WANTZ) Z( 1, 1 ) = ONE;
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

      ANRM = SLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      ISCALE = 0
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM > RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            slascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            slascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call SSBTRD to reduce symmetric band matrix to tridiagonal form.

      INDE = 1
      INDWRK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1
      ssbtrd(JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEDC.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         sstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         sgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N );
         slacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL SSCAL( N, ONE / SIGMA, W, 1 );

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of SSBEVD

      }
