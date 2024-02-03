      SUBROUTINE SLAQZ2( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ )
      IMPLICIT NONE

      // Arguments
      bool   , INTENT( IN ) :: ILQ, ILZ;
      int    , INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, NQ, NZ, QSTART, ZSTART, IHI;
      REAL :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )

      // Parameters
      REAL :: ZERO, ONE, HALF
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local variables
      REAL :: H( 2, 3 ), C1, S1, C2, S2, TEMP

      // External functions
      // EXTERNAL :: SLARTG, SROT

      IF( K+2 .EQ. IHI ) THEN
         // Shift is located on the edge of the matrix, remove it
         H = B( IHI-1:IHI, IHI-2:IHI )
         // Make H upper triangular
         CALL SLARTG( H( 1, 1 ), H( 2, 1 ), C1, S1, TEMP )
         H( 2, 1 ) = ZERO
         H( 1, 1 ) = TEMP
         CALL SROT( 2, H( 1, 2 ), 2, H( 2, 2 ), 2, C1, S1 )

         CALL SLARTG( H( 2, 3 ), H( 2, 2 ), C1, S1, TEMP )
         CALL SROT( 1, H( 1, 3 ), 1, H( 1, 2 ), 1, C1, S1 )
         CALL SLARTG( H( 1, 2 ), H( 1, 1 ), C2, S2, TEMP )

         CALL SROT( IHI-ISTARTM+1, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C1, S1 )          CALL SROT( IHI-ISTARTM+1, B( ISTARTM, IHI-1 ), 1, B( ISTARTM, IHI-2 ), 1, C2, S2 )
         B( IHI-1, IHI-2 ) = ZERO
         B( IHI, IHI-2 ) = ZERO
         CALL SROT( IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C1, S1 )          CALL SROT( IHI-ISTARTM+1, A( ISTARTM, IHI-1 ), 1, A( ISTARTM, IHI-2 ), 1, C2, S2 )
         IF ( ILZ ) THEN
            CALL SROT( NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C1, S1 )             CALL SROT( NZ, Z( 1, IHI-1-ZSTART+1 ), 1, Z( 1, IHI-2-ZSTART+1 ), 1, C2, S2 )
         END IF

         CALL SLARTG( A( IHI-1, IHI-2 ), A( IHI, IHI-2 ), C1, S1, TEMP )
         A( IHI-1, IHI-2 ) = TEMP
         A( IHI, IHI-2 ) = ZERO
         CALL SROT( ISTOPM-IHI+2, A( IHI-1, IHI-1 ), LDA, A( IHI, IHI-1 ), LDA, C1, S1 )          CALL SROT( ISTOPM-IHI+2, B( IHI-1, IHI-1 ), LDB, B( IHI, IHI-1 ), LDB, C1, S1 )
         IF ( ILQ ) THEN
            CALL SROT( NQ, Q( 1, IHI-1-QSTART+1 ), 1, Q( 1, IHI-QSTART+ 1 ), 1, C1, S1 )
         END IF

         CALL SLARTG( B( IHI, IHI ), B( IHI, IHI-1 ), C1, S1, TEMP )
         B( IHI, IHI ) = TEMP
         B( IHI, IHI-1 ) = ZERO
         CALL SROT( IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C1, S1 )          CALL SROT( IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C1, S1 )
         IF ( ILZ ) THEN
            CALL SROT( NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C1, S1 )
         END IF

      } else {

         // Normal operation, move bulge down

         H = B( K+1:K+2, K:K+2 )

         // Make H upper triangular

         CALL SLARTG( H( 1, 1 ), H( 2, 1 ), C1, S1, TEMP )
         H( 2, 1 ) = ZERO
         H( 1, 1 ) = TEMP
         CALL SROT( 2, H( 1, 2 ), 2, H( 2, 2 ), 2, C1, S1 )

         // Calculate Z1 and Z2

         CALL SLARTG( H( 2, 3 ), H( 2, 2 ), C1, S1, TEMP )
         CALL SROT( 1, H( 1, 3 ), 1, H( 1, 2 ), 1, C1, S1 )
         CALL SLARTG( H( 1, 2 ), H( 1, 1 ), C2, S2, TEMP )

         // Apply transformations from the right

         CALL SROT( K+3-ISTARTM+1, A( ISTARTM, K+2 ), 1, A( ISTARTM, K+1 ), 1, C1, S1 )          CALL SROT( K+3-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM, K ), 1, C2, S2 )          CALL SROT( K+2-ISTARTM+1, B( ISTARTM, K+2 ), 1, B( ISTARTM, K+1 ), 1, C1, S1 )          CALL SROT( K+2-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ), 1, C2, S2 )
         IF ( ILZ ) THEN
            CALL SROT( NZ, Z( 1, K+2-ZSTART+1 ), 1, Z( 1, K+1-ZSTART+ 1 ), 1, C1, S1 )             CALL SROT( NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ), 1, C2, S2 )
         END IF
         B( K+1, K ) = ZERO
         B( K+2, K ) = ZERO

         // Calculate Q1 and Q2

         CALL SLARTG( A( K+2, K ), A( K+3, K ), C1, S1, TEMP )
         A( K+2, K ) = TEMP
         A( K+3, K ) = ZERO
         CALL SLARTG( A( K+1, K ), A( K+2, K ), C2, S2, TEMP )
         A( K+1, K ) = TEMP
         A( K+2, K ) = ZERO

      // Apply transformations from the left

         CALL SROT( ISTOPM-K, A( K+2, K+1 ), LDA, A( K+3, K+1 ), LDA, C1, S1 )          CALL SROT( ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C2, S2 )

         CALL SROT( ISTOPM-K, B( K+2, K+1 ), LDB, B( K+3, K+1 ), LDB, C1, S1 )          CALL SROT( ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C2, S2 )
         IF ( ILQ ) THEN
            CALL SROT( NQ, Q( 1, K+2-QSTART+1 ), 1, Q( 1, K+3-QSTART+ 1 ), 1, C1, S1 )             CALL SROT( NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+ 1 ), 1, C2, S2 )
         END IF

      END IF

      // End of SLAQZ2

      END SUBROUTINE
