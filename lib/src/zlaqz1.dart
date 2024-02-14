      void zlaqz1(final int ILQ, final int ILZ, final int K, final int ISTARTM, final int ISTOPM, final int IHI, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int NQ, final int QSTART, final Matrix<double> Q_, final int LDQ, final int NZ, final int ZSTART, final int Z, final int LDZ,) {
  final A = A_.dim();
  final B = B_.dim();
  final Q = Q_.dim();
      // Arguments
      bool   , INTENT( IN ) :: ILQ, ILZ;
      int    , INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, NQ, NZ, QSTART, ZSTART, IHI;
      Complex :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );

      // Parameters
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local variables
      double           :: C;
      Complex :: S, TEMP;

      // External Functions
      // EXTERNAL :: ZLARTG, ZROT

      if ( K+1 == IHI ) {

         // Shift is located on the edge of the matrix, remove it

         zlartg(B( IHI, IHI ), B( IHI, IHI-1 ), C, S, TEMP );
         B[IHI][IHI] = TEMP;
         B[IHI][IHI-1] = CZERO;
         zrot(IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C, S );
         zrot(IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C, S );
         if ( ILZ ) {
            zrot(NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C, S );
         }

      } else {

         // Normal operation, move bulge down


         // Apply transformation from the right

         zlartg(B( K+1, K+1 ), B( K+1, K ), C, S, TEMP );
         B[K+1][K+1] = TEMP;
         B[K+1][K] = CZERO;
         zrot(K+2-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM, K ), 1, C, S );
         zrot(K-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ), 1, C, S );
         if ( ILZ ) {
            zrot(NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ), 1, C, S );
         }

         // Apply transformation from the left

         zlartg(A( K+1, K ), A( K+2, K ), C, S, TEMP );
         A[K+1][K] = TEMP;
         A[K+2][K] = CZERO;
         zrot(ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C, S );
         zrot(ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C, S );
         if ( ILQ ) {
            zrot(NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+ 1 ), 1, C, DCONJG( S ) );
         }

      }
      END SUBROUTINE;
