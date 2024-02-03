      SUBROUTINE ZPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      COMPLEX*16         Z( LDZ, * )
      // ..

*  ====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPTTRF, XERBLA, ZBDSQR, ZLASET
      // ..
      // .. Local Arrays ..
      COMPLEX*16         C( 1, 1 ), VT( 1, 1 )
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
         xerbla('ZPTEQR', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      if ( N.EQ.1 ) {
         if (ICOMPZ.GT.0) Z( 1, 1 ) = CONE;
         RETURN
      }
      if (ICOMPZ.EQ.2) CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ );

      // Call DPTTRF to factor the matrix.

      dpttrf(N, D, E, INFO );
      if (INFO.NE.0) RETURN;
      for (I = 1; I <= N; I++) { // 10
         D( I ) = SQRT( D( I ) )
      } // 10
      for (I = 1; I <= N - 1; I++) { // 20
         E( I ) = E( I )*D( I )
      } // 20

      // Call ZBDSQR to compute the singular values/vectors of the
      // bidiagonal factor.

      if ( ICOMPZ.GT.0 ) {
         NRU = N
      } else {
         NRU = 0
      }
      zbdsqr('Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO );

      // Square the singular values.

      if ( INFO.EQ.0 ) {
         for (I = 1; I <= N; I++) { // 30
            D( I ) = D( I )*D( I )
         } // 30
      } else {
         INFO = N + INFO
      }

      RETURN

      // End of ZPTEQR

      }
