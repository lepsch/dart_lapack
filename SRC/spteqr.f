      SUBROUTINE SPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBDSQR, SLASET, SPTTRF, XERBLA
      // ..
      // .. Local Arrays ..
      REAL               C( 1, 1 ), VT( 1, 1 )
      // ..
      // .. Local Scalars ..
      int                I, ICOMPZ, NRU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( LSAME( COMPZ, 'N' ) ) {
         ICOMPZ = 0
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ICOMPZ = 1
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ICOMPZ = 2
      } else {
         ICOMPZ = -1
      }
      if ( ICOMPZ.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('SPTEQR', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         IF( ICOMPZ.GT.0 ) Z( 1, 1 ) = ONE
         RETURN
      }
      IF( ICOMPZ.EQ.2 ) CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )

      // Call SPTTRF to factor the matrix.

      spttrf(N, D, E, INFO );
      IF( INFO.NE.0 ) RETURN
      for (I = 1; I <= N; I++) { // 10
         D( I ) = SQRT( D( I ) )
      } // 10
      DO 20 I = 1, N - 1
         E( I ) = E( I )*D( I )
      } // 20

      // Call SBDSQR to compute the singular values/vectors of the
      // bidiagonal factor.

      if ( ICOMPZ.GT.0 ) {
         NRU = N
      } else {
         NRU = 0
      }
      sbdsqr('Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO );

      // Square the singular values.

      if ( INFO.EQ.0 ) {
         for (I = 1; I <= N; I++) { // 30
            D( I ) = D( I )*D( I )
         } // 30
      } else {
         INFO = N + INFO
      }

      RETURN

      // End of SPTEQR

      }
