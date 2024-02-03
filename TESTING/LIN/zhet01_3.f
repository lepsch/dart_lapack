      SUBROUTINE ZHET01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID )
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
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J
      double             ANORM, EPS;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             ZLANHE, DLAMCH;
      EXTERNAL           LSAME, ZLANHE, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET, ZLAVHE_ROOK, ZSYCONVF_ROOK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIMAG, DBLE
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
*     a) Revert to multipliers of L
*
      CALL ZSYCONVF_ROOK( UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO )
*
*     1) Determine EPS and the norm of A.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      DO J = 1, N
         IF( DIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
      END DO
*
*     2) Initialize C to the identity matrix.
*
      CALL ZLASET( 'Full', N, N, CZERO, CONE, C, LDC )
*
*     3) Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).
*
      CALL ZLAVHE_ROOK( UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
*     4) Call ZLAVHE_RK again to multiply by U (or L ).
*
      CALL ZLAVHE_ROOK( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
*     5) Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            DO I = 1, J - 1
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
         END DO
      ELSE
         DO J = 1, N
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
            DO I = J + 1, N
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      END IF
*
*     6) Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID/DBLE( N ) )/ANORM ) / EPS
      END IF
*
*     b) Convert to factor of L (or U)
*
      CALL ZSYCONVF_ROOK( UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO )
*
      RETURN
*
*     End of ZHET01_3
*
      END
