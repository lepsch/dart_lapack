      SUBROUTINE CLAQZ1( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B,
     $                   LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ )
      IMPLICIT NONE
*
*     Arguments
      LOGICAL, INTENT( IN ) :: ILQ, ILZ
      INTEGER, INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM,
     $         NQ, NZ, QSTART, ZSTART, IHI
      COMPLEX :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
*
*     Parameters
      COMPLEX         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) )
      REAL :: ZERO, ONE, HALF
      PARAMETER( ZERO = 0.0, ONE = 1.0, HALF = 0.5 )
*
*     Local variables
      REAL :: C
      COMPLEX :: S, TEMP
*
*     External Functions
      EXTERNAL :: CLARTG, CROT
*
      IF( K+1 .EQ. IHI ) THEN
*
*        Shift is located on the edge of the matrix, remove it
*
         CALL CLARTG( B( IHI, IHI ), B( IHI, IHI-1 ), C, S, TEMP )
         B( IHI, IHI ) = TEMP
         B( IHI, IHI-1 ) = CZERO
         CALL CROT( IHI-ISTARTM, B( ISTARTM, IHI ), 1, B( ISTARTM,
     $              IHI-1 ), 1, C, S )
         CALL CROT( IHI-ISTARTM+1, A( ISTARTM, IHI ), 1, A( ISTARTM,
     $              IHI-1 ), 1, C, S )
         IF ( ILZ ) THEN
            CALL CROT( NZ, Z( 1, IHI-ZSTART+1 ), 1, Z( 1, IHI-1-ZSTART+
     $                 1 ), 1, C, S )
         END IF
*
      ELSE
*
*        Normal operation, move bulge down
*
*
*        Apply transformation from the right
*
         CALL CLARTG( B( K+1, K+1 ), B( K+1, K ), C, S, TEMP )
         B( K+1, K+1 ) = TEMP
         B( K+1, K ) = CZERO
         CALL CROT( K+2-ISTARTM+1, A( ISTARTM, K+1 ), 1, A( ISTARTM,
     $              K ), 1, C, S )
         CALL CROT( K-ISTARTM+1, B( ISTARTM, K+1 ), 1, B( ISTARTM, K ),
     $              1, C, S )
         IF ( ILZ ) THEN
            CALL CROT( NZ, Z( 1, K+1-ZSTART+1 ), 1, Z( 1, K-ZSTART+1 ),
     $                 1, C, S )
         END IF
*
*        Apply transformation from the left
*
         CALL CLARTG( A( K+1, K ), A( K+2, K ), C, S, TEMP )
         A( K+1, K ) = TEMP
         A( K+2, K ) = CZERO
         CALL CROT( ISTOPM-K, A( K+1, K+1 ), LDA, A( K+2, K+1 ), LDA, C,
     $              S )
         CALL CROT( ISTOPM-K, B( K+1, K+1 ), LDB, B( K+2, K+1 ), LDB, C,
     $              S )
         IF ( ILQ ) THEN
            CALL CROT( NQ, Q( 1, K+1-QSTART+1 ), 1, Q( 1, K+2-QSTART+
     $                 1 ), 1, C, CONJG( S ) )
         END IF
*
      END IF
*
*     End of CLAQZ1
*
      END SUBROUTINE