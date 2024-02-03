      SUBROUTINE DTPT01( UPLO, DIAG, N, AP, AINVP, RCOND, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             AINVP( * ), AP( * ), WORK( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               UNITD;
      int                J, JC;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANTP;
      // EXTERNAL LSAME, DLAMCH, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..
*
      // Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANTP( '1', UPLO, DIAG, N, AP, WORK )
      AINVNM = DLANTP( '1', UPLO, DIAG, N, AINVP, WORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
      // Compute A * AINV, overwriting AINV.
*
      UNITD = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         JC = 1
         DO 10 J = 1, N
            IF( UNITD ) AINVP( JC+J-1 ) = ONE
*
            // Form the j-th column of A*AINV
*
            CALL DTPMV( 'Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 )
*
            // Subtract 1 from the diagonal
*
            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE
            JC = JC + J
   10    CONTINUE
      ELSE
         JC = 1
         DO 20 J = 1, N
            IF( UNITD ) AINVP( JC ) = ONE
*
            // Form the j-th column of A*AINV
*
            CALL DTPMV( 'Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 )
*
            // Subtract 1 from the diagonal
*
            AINVP( JC ) = AINVP( JC ) - ONE
            JC = JC + N - J + 1
   20    CONTINUE
      END IF
*
      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = DLANTP( '1', UPLO, 'Non-unit', N, AINVP, WORK )
*
      RESID = ( ( RESID*RCOND ) / DBLE( N ) ) / EPS
*
      RETURN
*
      // End of DTPT01
*
      END
