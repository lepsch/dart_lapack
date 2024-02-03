      SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSP;
      // EXTERNAL LSAME, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DOPGTR, DSCAL, DSPTRD, DSTEQR, DSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( UPLO, 'U' ) || LSAME( UPLO, 'L' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDZ.LT.1 || ( WANTZ && LDZ.LT.N ) ) {
         INFO = -7
      }

      if ( INFO != 0 ) {
         xerbla('DSPEV ', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         W( 1 ) = AP( 1 )
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO && ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         dscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDE = 1
      INDTAU = INDE + N
      dsptrd(UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // DOPGTR to generate the orthogonal matrix, then call DSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         INDWRK = INDTAU + N
         dopgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
         dsteqr(JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDTAU ), INFO );
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

      RETURN

      // End of DSPEV

      }
