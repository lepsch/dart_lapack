      SUBROUTINE SSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                ISCALE, LIWMIN, LWMIN;
      REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TNRM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANST, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANST, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSTEDC, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK == -1 .OR. LIWORK == -1 )

      INFO = 0
      LIWMIN = 1
      LWMIN = 1
      if ( N.GT.1 && WANTZ ) {
         LWMIN = 1 + 4*N + N**2
         LIWMIN = 3 + 5*N
      }

      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.N ) ) {
         INFO = -6
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN && .NOT.LQUERY ) {
            INFO = -8
         } else if ( LIWORK.LT.LIWMIN && .NOT.LQUERY ) {
            INFO = -10
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSTEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         if (WANTZ) Z( 1, 1 ) = ONE;
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

      ISCALE = 0
      TNRM = SLANST( 'M', N, D, E )
      if ( TNRM.GT.ZERO && TNRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / TNRM
      } else if ( TNRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / TNRM
      }
      if ( ISCALE == 1 ) {
         sscal(N, SIGMA, D, 1 );
         sscal(N-1, SIGMA, E( 1 ), 1 );
      }

      // For eigenvalues only, call SSTERF.  For eigenvalues and
      // eigenvectors, call SSTEDC.

      if ( .NOT.WANTZ ) {
         ssterf(N, D, E, INFO );
      } else {
         sstedc('I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL SSCAL( N, ONE / SIGMA, D, 1 );

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SSTEVD

      }
