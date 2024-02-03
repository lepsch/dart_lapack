      SUBROUTINE ZTPT01( UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         AINVP( * ), AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UNITD
      int                J, JC
      DOUBLE PRECISION   AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANTP
      EXTERNAL           LSAME, DLAMCH, ZLANTP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANTP( '1', UPLO, DIAG, N, AP, RWORK )
      AINVNM = ZLANTP( '1', UPLO, DIAG, N, AINVP, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Compute A * AINV, overwriting AINV.
*
      UNITD = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         JC = 1
         DO 10 J = 1, N
            IF( UNITD ) AINVP( JC+J-1 ) = ONE
*
*           Form the j-th column of A*AINV.
*
            CALL ZTPMV( 'Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal to form A*AINV - I.
*
            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE
            JC = JC + J
   10    CONTINUE
      ELSE
         JC = 1
         DO 20 J = 1, N
            IF( UNITD ) AINVP( JC ) = ONE
*
*           Form the j-th column of A*AINV.
*
            CALL ZTPMV( 'Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal to form A*AINV - I.
*
            AINVP( JC ) = AINVP( JC ) - ONE
            JC = JC + N - J + 1
   20    CONTINUE
      END IF
*
*     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = ZLANTP( '1', UPLO, 'Non-unit', N, AINVP, RWORK )
*
      RESID = ( ( RESID*RCOND ) / DBLE( N ) ) / EPS
*
      RETURN
*
*     End of ZTPT01
*
      END
