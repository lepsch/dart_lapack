      SUBROUTINE DSBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N, LWORK;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, WANTZ, LQUERY;
      int                IINFO, IMAX, INDE, INDWRK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, DLANSB;
      // EXTERNAL LSAME, DLAMCH, DLANSB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DSCAL, DSTEQR, DSTERF, XERBLA, DSYTRD_SB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
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
         if ( N <= 1 ) {
            LWMIN = 1
            WORK( 1 ) = LWMIN
         } else {
            IB    = ILAENV2STAGE( 2, 'DSYTRD_SB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'DSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = N + LHTRD + LWTRD
            WORK( 1 )  = LWMIN
         }

         if (LWORK < LWMIN && .NOT.LQUERY) INFO = -11;
      }

      if ( INFO != 0 ) {
         xerbla('DSBEV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( LOWER ) {
            W( 1 ) = AB( 1, 1 )
         } else {
            W( 1 ) = AB( KD+1, 1 )
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = DLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
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
            dlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            dlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form.

      INDE    = 1
      INDHOUS = INDE + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      dsytrd_sb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dsteqr(JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = LWMIN

      RETURN

      // End of DSBEV_2STAGE

      }
