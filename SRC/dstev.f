      SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               WANTZ;
      int                IMAX, ISCALE;
      double             BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TNRM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTEQR, DSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -6
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DSTEV ', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         IF( WANTZ ) Z( 1, 1 ) = ONE
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
         CALL DSCAL( N, SIGMA, D, 1 )
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
      }

      // For eigenvalues only, call DSTERF.  For eigenvalues and
      // eigenvectors, call DSTEQR.

      if ( .NOT.WANTZ ) {
         CALL DSTERF( N, D, E, INFO )
      } else {
         CALL DSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         CALL DSCAL( IMAX, ONE / SIGMA, D, 1 )
      }

      RETURN

      // End of DSTEV

      }
