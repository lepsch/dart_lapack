      SUBROUTINE DSYT01( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
*     ..
*     .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J;
      double             ANORM, EPS;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSY;
      EXTERNAL           LSAME, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET, DLAVSY
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE
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
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )
*
*     Initialize C to the identity matrix.
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, C, LDC )
*
*     Call DLAVSY to form the product D * U' (or D * L' ).
*
      CALL DLAVSY( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
*     Call DLAVSY again to multiply by U (or L ).
*
      CALL DLAVSY( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
*     Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J
               C( I, J ) = C( I, J ) - A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = J, N
               C( I, J ) = C( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = DLANSY( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of DSYT01
*
      END
