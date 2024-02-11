      void zpteqr(final int COMPZ, final int N, final int D, final int E, final Matrix<double> Z, final int LDZ, final Array<double> _WORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             COMPZ;
      int                INFO, LDZ, N;
      double             D( * ), E( * ), WORK( * );
      Complex         Z( LDZ, * );
      // ..

// ====================================================================

      // .. Parameters ..
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPTTRF, XERBLA, ZBDSQR, ZLASET
      Complex         C( 1, 1 ), VT( 1, 1 );
      int                I, ICOMPZ, NRU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters.

      INFO = 0;

      if ( lsame( COMPZ, 'N' ) ) {
         ICOMPZ = 0;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ICOMPZ = 2;
      } else {
         ICOMPZ = -1;
      }
      if ( ICOMPZ < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ( LDZ < 1 ) || ( ICOMPZ > 0 && LDZ < max( 1, N ) ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('ZPTEQR', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if (ICOMPZ > 0) Z( 1, 1 ) = CONE;
         return;
      }
      if (ICOMPZ == 2) zlaset( 'Full', N, N, CZERO, CONE, Z, LDZ );

      // Call DPTTRF to factor the matrix.

      dpttrf(N, D, E, INFO );
      if (INFO != 0) return;
      for (I = 1; I <= N; I++) { // 10
         D[I] = sqrt( D( I ) );
      } // 10
      for (I = 1; I <= N - 1; I++) { // 20
         E[I] = E( I )*D( I );
      } // 20

      // Call ZBDSQR to compute the singular values/vectors of the
      // bidiagonal factor.

      if ( ICOMPZ > 0 ) {
         NRU = N;
      } else {
         NRU = 0;
      }
      zbdsqr('Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO );

      // Square the singular values.

      if ( INFO == 0 ) {
         for (I = 1; I <= N; I++) { // 30
            D[I] = D( I )*D( I );
         } // 30
      } else {
         INFO = N + INFO;
      }

      }
