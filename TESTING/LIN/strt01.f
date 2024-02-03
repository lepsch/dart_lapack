      SUBROUTINE STRT01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                LDA, LDAINV, N
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AINV( LDAINV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                J
      REAL               AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANTR
      EXTERNAL           LSAME, SLAMCH, SLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           STRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0
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
      ANORM = SLANTR( '1', UPLO, DIAG, N, N, A, LDA, WORK )
      AINVNM = SLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, WORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Set the diagonal of AINV to 1 if AINV has unit diagonal.
*
      IF( LSAME( DIAG, 'U' ) ) THEN
         DO 10 J = 1, N
            AINV( J, J ) = ONE
   10    CONTINUE
      END IF
*
*     Compute A * AINV, overwriting AINV.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            CALL STRMV( 'Upper', 'No transpose', DIAG, J, A, LDA, AINV( 1, J ), 1 )
   20    CONTINUE
      ELSE
         DO 30 J = 1, N
            CALL STRMV( 'Lower', 'No transpose', DIAG, N-J+1, A( J, J ), LDA, AINV( J, J ), 1 )
   30    CONTINUE
      END IF
*
*     Subtract 1 from each diagonal element to form A*AINV - I.
*
      DO 40 J = 1, N
         AINV( J, J ) = AINV( J, J ) - ONE
   40 CONTINUE
*
*     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = SLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, WORK )
*
      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS
*
      RETURN
*
*     End of STRT01
*
      END
