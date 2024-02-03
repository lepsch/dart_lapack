      SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
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
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, INDWRK, ISCALE, LIOPT, LIWMIN, LLRWK, LLWORK, LLWRK2, LOPT, LROPT, LRWMIN, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANHE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZHETRD, ZLACPY, ZLASCL, ZSTEDC, ZUNMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

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
         if ( N.LE.1 ) {
            LWMIN = 1
            LRWMIN = 1
            LIWMIN = 1
            LOPT = LWMIN
            LROPT = LRWMIN
            LIOPT = LIWMIN
         } else {
            if ( WANTZ ) {
               LWMIN = 2*N + N*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            } else {
               LWMIN = N + 1
               LRWMIN = N
               LIWMIN = 1
            }
            LOPT = MAX( LWMIN, N + N*ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
            LROPT = LRWMIN
            LIOPT = LIWMIN
         }
         WORK( 1 ) = LOPT
         RWORK( 1 ) = LROPT
         IWORK( 1 ) = LIOPT

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -8
         } else if ( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) {
            INFO = -10
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHEEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         W( 1 ) = DBLE( A( 1, 1 ) )
         IF( WANTZ ) A( 1, 1 ) = CONE
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
      IF( ISCALE.EQ.1 ) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )

      // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      INDRWK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWORK = LWORK - INDWRK + 1
      LLWRK2 = LWORK - INDWK2 + 1
      LLRWK = LRWORK - INDRWK + 1
      zhetrd(UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
      // tridiagonal matrix, then call ZUNMTR to multiply it to the
      // Householder transformations represented as Householder vectors in
      // A.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         zunmtr('L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO );
         zlacpy('A', N, N, WORK( INDWRK ), N, A, LDA );
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

      WORK( 1 ) = LOPT
      RWORK( 1 ) = LROPT
      IWORK( 1 ) = LIOPT

      RETURN

      // End of ZHEEVD

      }
