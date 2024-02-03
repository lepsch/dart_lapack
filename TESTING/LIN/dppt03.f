      SUBROUTINE DPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                LDWORK, N
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), AINV( * ), RWORK( * ), WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J, JJ
      DOUBLE PRECISION   AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSP
      EXTERNAL           LSAME, DLAMCH, DLANGE, DLANSP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSPMV
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
      ANORM = DLANSP( '1', UPLO, N, A, RWORK )
      AINVNM = DLANSP( '1', UPLO, N, AINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.EQ.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     UPLO = 'U':
*     Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
*     expand it to a full matrix, then multiply by A one column at a
*     time, moving the result one column to the left.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Copy AINV
*
         JJ = 1
         DO 10 J = 1, N - 1
            CALL DCOPY( J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 )
            CALL DCOPY( J-1, AINV( JJ ), 1, WORK( J, 2 ), LDWORK )
            JJ = JJ + J
   10    CONTINUE
         JJ = ( ( N-1 )*N ) / 2 + 1
         CALL DCOPY( N-1, AINV( JJ ), 1, WORK( N, 2 ), LDWORK )
*
*        Multiply by A
*
         DO 20 J = 1, N - 1
            CALL DSPMV( 'Upper', N, -ONE, A, WORK( 1, J+1 ), 1, ZERO, WORK( 1, J ), 1 )
   20    CONTINUE
         CALL DSPMV( 'Upper', N, -ONE, A, AINV( JJ ), 1, ZERO, WORK( 1, N ), 1 )
*
*     UPLO = 'L':
*     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
*     and multiply by A, moving each column to the right.
*
      ELSE
*
*        Copy AINV
*
         CALL DCOPY( N-1, AINV( 2 ), 1, WORK( 1, 1 ), LDWORK )
         JJ = N + 1
         DO 30 J = 2, N
            CALL DCOPY( N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 )
            CALL DCOPY( N-J, AINV( JJ+1 ), 1, WORK( J, J ), LDWORK )
            JJ = JJ + N - J + 1
   30    CONTINUE
*
*        Multiply by A
*
         DO 40 J = N, 2, -1
            CALL DSPMV( 'Lower', N, -ONE, A, WORK( 1, J-1 ), 1, ZERO, WORK( 1, J ), 1 )
   40    CONTINUE
         CALL DSPMV( 'Lower', N, -ONE, A, AINV( 1 ), 1, ZERO, WORK( 1, 1 ), 1 )
*
      END IF
*
*     Add the identity matrix to WORK .
*
      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + ONE
   50 CONTINUE
*
*     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = DLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )
*
      RETURN
*
*     End of DPPT03
*
      END
