      void claqz1(ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ ) {
      // IMPLICIT NONE

      // Arguments
      bool   , INTENT( IN ) :: ILQ, ILZ;
      int    , INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, NQ, NZ, QSTART, ZSTART, IHI;
      COMPLEX :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );

      // Parameters
      COMPLEX         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      REAL :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local variables
      REAL :: C;
      COMPLEX :: S, TEMP;

      // External Functions
      // EXTERNAL :: CLARTG, CROT

      if ( K+1 == IHI ) {

         // Shift is located on the edge of the matrix, remove it

         clartg(B( IHI, IHI ), B( IHI, IHI-1 ), C, S, TEMP );
         B( IHI, IHI ) = TEMP;
         B( IHI, IHI-1 ) = CZERO;
         crot(IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C, S );
         crot(IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C, S );
         if ( ILZ ) {
            crot(NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C, S );
         }

      } else {

         // Normal operation, move bulge down


         // Apply transformation from the right

         clartg(B( K+1, K+1 ), B( K+1, K ), C, S, TEMP );
         B( K+1, K+1 ) = TEMP;
         B( K+1, K ) = CZERO;
         crot(K+2-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM, K ), 1, C, S );
         crot(K-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ), 1, C, S );
         if ( ILZ ) {
            crot(NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ), 1, C, S );
         }

         // Apply transformation from the left

         clartg(A( K+1, K ), A( K+2, K ), C, S, TEMP );
         A( K+1, K ) = TEMP;
         A( K+2, K ) = CZERO;
         crot(ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C, S );
         crot(ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C, S );
         if ( ILQ ) {
            crot(NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+ 1 ), 1, C, CONJG( S ) );
         }

      }
      END SUBROUTINE;
