      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), W( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWKOPT, NB;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANHE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZHETRD, ZLASCL, ZSTEQR, ZUNGTR
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
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }

      if ( INFO.EQ.0 ) {
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+1 )*N )
         WORK( 1 ) = LWKOPT

         IF( LWORK.LT.MAX( 1, 2*N-1 ) .AND. .NOT.LQUERY ) INFO = -8
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHEEV ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         RETURN
      }

      if ( N.EQ.1 ) {
         W( 1 ) = DBLE( A( 1, 1 ) )
         WORK( 1 ) = 1
         if (WANTZ) A( 1, 1 ) = CONE;
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

      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if (ISCALE.EQ.1) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO );

      // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      zhetrd(UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // ZUNGTR to generate the unitary matrix, then call ZSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zungtr(UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
         INDWRK = INDE + N
         zsteqr(JOBZ, N, W, RWORK( INDE ), A, LDA, RWORK( INDWRK ), INFO );
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

      // Set WORK(1) to optimal complex workspace size.

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZHEEV

      }
