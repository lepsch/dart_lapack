      void dstevd(JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTZ;
      int                ISCALE, LIWMIN, LWMIN;
      double             BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TNRM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTEDC, DSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      LIWMIN = 1;
      LWMIN = 1;
      if ( N > 1 && WANTZ ) {
         LWMIN = 1 + 4*N + N**2;
         LIWMIN = 3 + 5*N;
      }

      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -6;
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = LWMIN;
         IWORK( 1 ) = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -10;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSTEVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if (WANTZ) Z( 1, 1 ) = ONE;
         return;
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = sqrt( BIGNUM );

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0;
      TNRM = DLANST( 'M', N, D, E );
      if ( TNRM > ZERO && TNRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / TNRM;
      } else if ( TNRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / TNRM;
      }
      if ( ISCALE == 1 ) {
         dscal(N, SIGMA, D, 1 );
         dscal(N-1, SIGMA, E( 1 ), 1 );
      }

      // For eigenvalues only, call DSTERF.  For eigenvalues and
      // eigenvectors, call DSTEDC.

      if ( !WANTZ ) {
         dsterf(N, D, E, INFO );
      } else {
         dstedc('I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE == 1) CALL DSCAL( N, ONE / SIGMA, D, 1 );

      WORK( 1 ) = LWMIN;
      IWORK( 1 ) = LIWMIN;

      return;

      // End of DSTEVD

      }
