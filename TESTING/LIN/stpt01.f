      SUBROUTINE STPT01( UPLO, DIAG, N, AP, AINVP, RCOND, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            N
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               AINVP( * ), AP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UNITD
      INTEGER            J, JC
      REAL               AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANTP
      EXTERNAL           LSAME, SLAMCH, SLANTP
*     ..
*     .. External Subroutines ..
      EXTERNAL           STPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
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
      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANTP( '1', UPLO, DIAG, N, AP, WORK )
      AINVNM = SLANTP( '1', UPLO, DIAG, N, AINVP, WORK )
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
*           Form the j-th column of A*AINV
*
            CALL STPMV( 'Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal
*
            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE
            JC = JC + J
   10    CONTINUE
      ELSE
         JC = 1
         DO 20 J = 1, N
            IF( UNITD ) AINVP( JC ) = ONE
*
*           Form the j-th column of A*AINV
*
            CALL STPMV( 'Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal
*
            AINVP( JC ) = AINVP( JC ) - ONE
            JC = JC + N - J + 1
   20    CONTINUE
      END IF
*
*     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = SLANTP( '1', UPLO, 'Non-unit', N, AINVP, WORK )
*
      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS
*
      RETURN
*
*     End of STPT01
*
      END
