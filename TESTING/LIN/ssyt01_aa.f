      SUBROUTINE SSYT01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      REAL               RESID
*     ..
*     .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSY
      // EXTERNAL LSAME, SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      // EXTERNAL SLASET, SLAVSY, SSWAP, STRMM, SLACPY
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC REAL
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
      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
*
*     Initialize C to the tridiagonal matrix T.
*
      CALL SLASET( 'Full', N, N, ZERO, ZERO, C, LDC )
      CALL SLACPY( 'F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 )
      IF( N.GT.1 ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL SLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL SLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         ELSE
            CALL SLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL SLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         ENDIF
*
*        Call STRMM to form the product U' * D (or L * D ).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL STRMM( 'Left', UPLO, 'Transpose', 'Unit', N-1, N, ONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC )
         ELSE
            CALL STRMM( 'Left', UPLO, 'No transpose', 'Unit', N-1, N, ONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC )
         END IF
*
*        Call STRMM again to multiply by U (or L ).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL STRMM( 'Right', UPLO, 'No transpose', 'Unit', N, N-1, ONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC )
         ELSE
            CALL STRMM( 'Right', UPLO, 'Transpose', 'Unit', N, N-1, ONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC )
         END IF
      ENDIF
*
*     Apply symmetric pivots
*
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL SSWAP( N, C( J, 1 ), LDC, C( I, 1 ), LDC )
      END DO
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL SSWAP( N, C( 1, J ), 1, C( 1, I ), 1 )
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
      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of SSYT01_AA
*
      END
