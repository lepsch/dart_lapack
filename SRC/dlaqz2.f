      SUBROUTINE DLAQZ2( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ )
      IMPLICIT NONE

      // Arguments
      bool   , INTENT( IN ) :: ILQ, ILZ;
      int    , INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, NQ, NZ, QSTART, ZSTART, IHI;
      double           :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );

      // Parameters
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 ;

      // Local variables
      double           :: H( 2, 3 ), C1, S1, C2, S2, TEMP;

      // External functions
      // EXTERNAL :: DLARTG, DROT

      if ( K+2 .EQ. IHI ) {
         // Shift is located on the edge of the matrix, remove it
         H = B( IHI-1:IHI, IHI-2:IHI )
         // Make H upper triangular
         dlartg(H( 1, 1 ), H( 2, 1 ), C1, S1, TEMP );
         H( 2, 1 ) = ZERO
         H( 1, 1 ) = TEMP
         drot(2, H( 1, 2 ), 2, H( 2, 2 ), 2, C1, S1 );

         dlartg(H( 2, 3 ), H( 2, 2 ), C1, S1, TEMP );
         drot(1, H( 1, 3 ), 1, H( 1, 2 ), 1, C1, S1 );
         dlartg(H( 1, 2 ), H( 1, 1 ), C2, S2, TEMP );

         drot(IHI-ISTARTM+1, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C1, S1 );
         drot(IHI-ISTARTM+1, B( ISTARTM, IHI-1 ), 1, B( ISTARTM, IHI-2 ), 1, C2, S2 );
         B( IHI-1, IHI-2 ) = ZERO
         B( IHI, IHI-2 ) = ZERO
         drot(IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C1, S1 );
         drot(IHI-ISTARTM+1, A( ISTARTM, IHI-1 ), 1, A( ISTARTM, IHI-2 ), 1, C2, S2 );
         if ( ILZ ) {
            drot(NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C1, S1 );
            drot(NZ, Z( 1, IHI-1-ZSTART+1 ), 1, Z( 1, IHI-2-ZSTART+1 ), 1, C2, S2 );
         }

         dlartg(A( IHI-1, IHI-2 ), A( IHI, IHI-2 ), C1, S1, TEMP );
         A( IHI-1, IHI-2 ) = TEMP
         A( IHI, IHI-2 ) = ZERO
         drot(ISTOPM-IHI+2, A( IHI-1, IHI-1 ), LDA, A( IHI, IHI-1 ), LDA, C1, S1 );
         drot(ISTOPM-IHI+2, B( IHI-1, IHI-1 ), LDB, B( IHI, IHI-1 ), LDB, C1, S1 );
         if ( ILQ ) {
            drot(NQ, Q( 1, IHI-1-QSTART+1 ), 1, Q( 1, IHI-QSTART+ 1 ), 1, C1, S1 );
         }

         dlartg(B( IHI, IHI ), B( IHI, IHI-1 ), C1, S1, TEMP );
         B( IHI, IHI ) = TEMP
         B( IHI, IHI-1 ) = ZERO
         drot(IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C1, S1 );
         drot(IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C1, S1 );
         if ( ILZ ) {
            drot(NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C1, S1 );
         }

      } else {

         // Normal operation, move bulge down

         H = B( K+1:K+2, K:K+2 )

         // Make H upper triangular

         dlartg(H( 1, 1 ), H( 2, 1 ), C1, S1, TEMP );
         H( 2, 1 ) = ZERO
         H( 1, 1 ) = TEMP
         drot(2, H( 1, 2 ), 2, H( 2, 2 ), 2, C1, S1 );

         // Calculate Z1 and Z2

         dlartg(H( 2, 3 ), H( 2, 2 ), C1, S1, TEMP );
         drot(1, H( 1, 3 ), 1, H( 1, 2 ), 1, C1, S1 );
         dlartg(H( 1, 2 ), H( 1, 1 ), C2, S2, TEMP );

         // Apply transformations from the right

         drot(K+3-ISTARTM+1, A( ISTARTM, K+2 ), 1, A( ISTARTM, K+1 ), 1, C1, S1 );
         drot(K+3-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM, K ), 1, C2, S2 );
         drot(K+2-ISTARTM+1, B( ISTARTM, K+2 ), 1, B( ISTARTM, K+1 ), 1, C1, S1 );
         drot(K+2-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ), 1, C2, S2 );
         if ( ILZ ) {
            drot(NZ, Z( 1, K+2-ZSTART+1 ), 1, Z( 1, K+1-ZSTART+ 1 ), 1, C1, S1 );
            drot(NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ), 1, C2, S2 );
         }
         B( K+1, K ) = ZERO
         B( K+2, K ) = ZERO

         // Calculate Q1 and Q2

         dlartg(A( K+2, K ), A( K+3, K ), C1, S1, TEMP );
         A( K+2, K ) = TEMP
         A( K+3, K ) = ZERO
         dlartg(A( K+1, K ), A( K+2, K ), C2, S2, TEMP );
         A( K+1, K ) = TEMP
         A( K+2, K ) = ZERO

         // Apply transformations from the left

         drot(ISTOPM-K, A( K+2, K+1 ), LDA, A( K+3, K+1 ), LDA, C1, S1 );
         drot(ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C2, S2 );

         drot(ISTOPM-K, B( K+2, K+1 ), LDB, B( K+3, K+1 ), LDB, C1, S1 );
         drot(ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C2, S2 );
         if ( ILQ ) {
            drot(NQ, Q( 1, K+2-QSTART+1 ), 1, Q( 1, K+3-QSTART+ 1 ), 1, C1, S1 );
            drot(NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+ 1 ), 1, C2, S2 );
         }

      }

      // End of DLAQZ2

      END SUBROUTINE
