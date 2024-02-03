      SUBROUTINE CHEEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, CLANHE, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, CLANHE, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTERF, XERBLA, CLASCL, CSTEQR, CUNGTR, CHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, MAX, SQRT
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
         KD    = ILAENV2STAGE( 1, 'CHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'CHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = N + LHTRD + LWTRD
         WORK( 1 )  = SROUNDUP_LWORK(LWMIN)

         if (LWORK.LT.LWMIN .AND. .NOT.LQUERY) INFO = -8;
      }

      if ( INFO.NE.0 ) {
         xerbla('CHEEV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         RETURN
      }

      if ( N.EQ.1 ) {
         W( 1 ) = REAL( A( 1, 1 ) )
         WORK( 1 ) = 1
         if (WANTZ) A( 1, 1 ) = CONE;
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

      ANRM = CLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if (ISCALE.EQ.1) CALL CLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDE    = 1
      INDTAU  = 1
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      chetrd_2stage(JOBZ, UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // CUNGTR to generate the unitary matrix, then call CSTEQR.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cungtr(UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
         INDWRK = INDE + N
         csteqr(JOBZ, N, W, RWORK( INDE ), A, LDA, RWORK( INDWRK ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         sscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal complex workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of CHEEV_2STAGE

      }
