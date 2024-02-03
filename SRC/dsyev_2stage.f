      SUBROUTINE DSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), W( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, XERBLA, DSYTRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }

      if ( INFO.EQ.0 ) {
         KD    = ILAENV2STAGE( 1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = 2*N + LHTRD + LWTRD
         WORK( 1 )  = LWMIN

         if (LWORK.LT.LWMIN .AND. .NOT.LQUERY) INFO = -8;
      }

      if ( INFO.NE.0 ) {
         xerbla('DSYEV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         RETURN
      }

      if ( N.EQ.1 ) {
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         if (WANTZ) A( 1, 1 ) = ONE;
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

      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if (ISCALE.EQ.1) CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.

      INDE    = 1
      INDTAU  = INDE + N
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      dsytrd_2stage(JOBZ, UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // DORGTR to generate the orthogonal matrix, then call DSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         // Not available in this release, and argument checking should not
         // let it getting here
         RETURN
         dorgtr(UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )          CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = LWMIN

      RETURN

      // End of DSYEV_2STAGE

      }
