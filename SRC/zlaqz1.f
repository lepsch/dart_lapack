      SUBROUTINE ZLAQZ1( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ )
      IMPLICIT NONE
*
*     Arguments
      LOGICAL, INTENT( IN ) :: ILQ, ILZ
      int    , INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, NQ, NZ, QSTART, ZSTART, IHI
      COMPLEX*16 :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
*
*     Parameters
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION :: ZERO, ONE, HALF
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
*
*     Local variables
      DOUBLE PRECISION :: C
      COMPLEX*16 :: S, TEMP
*
*     External Functions
      EXTERNAL :: ZLARTG, ZROT
*
      IF( K+1 .EQ. IHI ) THEN
*
*        Shift is located on the edge of the matrix, remove it
*
         CALL ZLARTG( B( IHI, IHI ), B( IHI, IHI-1 ), C, S, TEMP )
         B( IHI, IHI ) = TEMP
         B( IHI, IHI-1 ) = CZERO
         CALL ZROT( IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM, IHI-1 ), 1, C, S )          CALL ZROT( IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM, IHI-1 ), 1, C, S )
         IF ( ILZ ) THEN
            CALL ZROT( NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+ 1 ), 1, C, S )
         END IF
*
      ELSE
*
*        Normal operation, move bulge down
*
*
*        Apply transformation from the right
*
         CALL ZLARTG( B( K+1, K+1 ), B( K+1, K ), C, S, TEMP )
         B( K+1, K+1 ) = TEMP
         B( K+1, K ) = CZERO
         CALL ZROT( K+2-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM, K ), 1, C, S )          CALL ZROT( K-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ), 1, C, S )
         IF ( ILZ ) THEN
            CALL ZROT( NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ), 1, C, S )
         END IF
*
*        Apply transformation from the left
*
         CALL ZLARTG( A( K+1, K ), A( K+2, K ), C, S, TEMP )
         A( K+1, K ) = TEMP
         A( K+2, K ) = CZERO
         CALL ZROT( ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C, S )          CALL ZROT( ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C, S )
         IF ( ILQ ) THEN
            CALL ZROT( NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+ 1 ), 1, C, DCONJG( S ) )
         END IF
*
      END IF
*
*     End of ZLAQZ1
*
      END SUBROUTINE
