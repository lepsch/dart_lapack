      SUBROUTINE ZHEEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )

      IMPLICIT NONE

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
      int                IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS       double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, ZLANHE;
      // EXTERNAL LSAME, DLAMCH, ZLANHE, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZLASCL, ZSTEQR, ZUNGTR, ZHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF

      IF( INFO.EQ.0 ) THEN
         KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = N + LHTRD + LWTRD
         WORK( 1 )  = LWMIN

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -8
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEV_2STAGE ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         RETURN
      END IF

      IF( N.EQ.1 ) THEN
         W( 1 ) = DBLE( A( 1, 1 ) )
         WORK( 1 ) = 1
         IF( WANTZ ) A( 1, 1 ) = CONE
         RETURN
      END IF

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )

      // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDE    = 1
      INDTAU  = 1
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      CALL ZHETRD_2STAGE( JOBZ, UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // ZUNGTR to generate the unitary matrix, then call ZSTEQR.

      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZUNGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
         INDWRK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), A, LDA, RWORK( INDWRK ), INFO )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF

      // Set WORK(1) to optimal complex workspace size.

      WORK( 1 ) = LWMIN

      RETURN

      // End of ZHEEV_2STAGE

      }
