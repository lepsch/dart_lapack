      SUBROUTINE CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            AP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK, ISCALE, LIWMIN, LLRWK, LLWRK, LRWMIN, LWMIN;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, CLANHP, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPTRD, CSSCAL, CSTEDC, CUPMTR, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK == -1 .OR. LRWORK == -1 .OR. LIWORK == -1 )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.N ) ) {
         INFO = -7
      }

      if ( INFO == 0 ) {
         if ( N.LE.1 ) {
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
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN && .NOT.LQUERY ) {
            INFO = -9
         } else if ( LRWORK.LT.LRWMIN && .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK.LT.LIWMIN && .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHPEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         W( 1 ) = REAL( AP( 1 ) )
         if (WANTZ) Z( 1, 1 ) = CONE;
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

      ANRM = CLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO && ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         csscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
      }

      // Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form.

      INDE = 1
      INDTAU = 1
      INDRWK = INDE + N
      INDWRK = INDTAU + N
      LLWRK = LWORK - INDWRK + 1
      LLRWK = LRWORK - INDRWK + 1
      chptrd(UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ), IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, first call
      // CUPGTR to generate the orthogonal matrix, then call CSTEDC.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cstedc('I', N, W, RWORK( INDE ), Z, LDZ, WORK( INDWRK ), LLWRK, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         cupmtr('L', UPLO, 'N', N, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
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

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of CHPEVD

      }
