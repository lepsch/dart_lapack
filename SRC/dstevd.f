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
      IF( N.GT.1 .AND. WANTZ ) THEN
         LWMIN = 1 + 4*N + N**2
         LIWMIN = 3 + 5*N
      END IF

      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -6
      END IF

      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      END IF

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
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, SIGMA, D, 1 )
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
      END IF

      // For eigenvalues only, call DSTERF.  For eigenvalues and
      // eigenvectors, call DSTEDC.

      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, D, E, INFO )
      } else {
         CALL DSTEDC( 'I', N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

      IF( ISCALE.EQ.1 ) CALL DSCAL( N, ONE / SIGMA, D, 1 )

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DSTEVD

      }
