      SUBROUTINE ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK, ISCALE, LIWMIN, LLRWK, LLWRK, LRWMIN, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHP;
      // EXTERNAL LSAME, DLAMCH, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZDSCAL, ZHPTRD, ZSTEDC, ZUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( UPLO, 'L' ) || LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -7
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         } else {
            if ( WANTZ ) {
               LWMIN = 2*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            } else {
               LWMIN = N
               LRWMIN = N
               LIWMIN = 1
            }
         }
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -9
         } else if ( LRWORK < LRWMIN && .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHPEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         W( 1 ) = DBLE( AP( 1 ) )
         if (WANTZ) Z( 1, 1 ) = CONE;
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

      ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM > RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         zdscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.

      INDE = 1
      INDTAU = 1
      INDRWK = INDE + N
      INDWRK = INDTAU + N
      LLWRK = LWORK - INDWRK + 1
      LLRWK = LRWORK - INDRWK + 1
      zhptrd(UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // ZUPGTR to generate the orthogonal matrix, then call ZSTEDC.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWRK, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         zupmtr('L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
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

      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of ZHPEVD

      }
