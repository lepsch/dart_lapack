      SUBROUTINE ZSYT01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N
      double             RESID;
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
      double             RWORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16          CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE  = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, J
      double             ANORM, EPS;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANSY;
      EXTERNAL           LSAME, DLAMCH, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET, ZLAVSY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Determine EPS and the norm of A.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK )
*
*     Initialize C to the tridiagonal matrix T.
*
      CALL ZLASET( 'Full', N, N, CZERO, CZERO, C, LDC )
      CALL ZLACPY( 'F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 )
      IF( N.GT.1 ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL ZLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL ZLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         ELSE
            CALL ZLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL ZLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         ENDIF
*
*        Call ZTRMM to form the product U' * D (or L * D ).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL ZTRMM( 'Left', UPLO, 'Transpose', 'Unit', N-1, N, CONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC )
         ELSE
            CALL ZTRMM( 'Left', UPLO, 'No transpose', 'Unit', N-1, N, CONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC )
         END IF
*
*        Call ZTRMM again to multiply by U (or L ).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL ZTRMM( 'Right', UPLO, 'No transpose', 'Unit', N, N-1, CONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC )
         ELSE
            CALL ZTRMM( 'Right', UPLO, 'Transpose', 'Unit', N, N-1, CONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC )
         END IF
      ENDIF
*
*     Apply symmetric pivots
*
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL ZSWAP( N, C( J, 1 ), LDC, C( I, 1 ), LDC )
      END DO
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL ZSWAP( N, C( 1, J ), 1, C( 1, I ), 1 )
      END DO
*
*
*     Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            DO I = 1, J
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      ELSE
         DO J = 1, N
            DO I = J, N
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      END IF
*
*     Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = ZLANSY( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of ZSYT01_AA
*
      END
