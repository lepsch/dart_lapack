      SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
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

      WANTZ = LSAME( JOBZ, 'V' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      LIWMIN = 1
      LWMIN = 1
      if ( N.GT.1 .AND. WANTZ ) {
         LWMIN = 1 + 4*N + N**2
         LIWMIN = 3 + 5*N
      }

      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -6
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -8
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -10
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DSTEVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      if ( N.EQ.1 ) {
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

      ISCALE = 0
      TNRM = DLANST( 'M', N, D, E )
      if ( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / TNRM
      } else if ( TNRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / TNRM
      }
      if ( ISCALE.EQ.1 ) {
         dscal(N, SIGMA, D, 1 );
         dscal(N-1, SIGMA, E( 1 ), 1 );
      }

      // For eigenvalues only, call DSTERF.  For eigenvalues and
      // eigenvectors, call DSTEDC.

      if ( .NOT.WANTZ ) {
         dsterf(N, D, E, INFO );
      } else {
         dstedc('I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE.EQ.1) CALL DSCAL( N, ONE / SIGMA, D, 1 );

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DSTEVD

      }
