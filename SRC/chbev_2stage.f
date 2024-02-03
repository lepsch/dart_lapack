      SUBROUTINE CHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N, LWORK;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, WANTZ, LQUERY;
      int                IINFO, IMAX, INDE, INDWRK, INDRWK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, CLANHB, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, CLANHB, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTERF, XERBLA, CLASCL, CSTEQR, CHETRD_2STAGE, CHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK == -1 )

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
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
         if ( N.LE.1 ) {
            LWMIN = 1
            WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         } else {
            IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = LHTRD + LWTRD
            WORK( 1 )  = SROUNDUP_LWORK(LWMIN)
         }

         if (LWORK < LWMIN && .NOT.LQUERY) INFO = -11;
      }

      if ( INFO != 0 ) {
         xerbla('CHBEV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( LOWER ) {
            W( 1 ) = REAL( AB( 1, 1 ) )
         } else {
            W( 1 ) = REAL( AB( KD+1, 1 ) )
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
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
            clascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            clascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDE    = 1
      INDHOUS = 1
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      chetrd_hb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         INDRWK = INDE + N
         csteqr(JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDRWK ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         sscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of CHBEV_2STAGE

      }
