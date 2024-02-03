      SUBROUTINE SPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAINV, LDWORK, N;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               AINVNM, ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL LSAME, SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RCOND = ONE
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
      AINVNM = SLANSY( '1', UPLO, N, AINV, LDAINV, RWORK )
      if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      }
      RCOND = ( ONE / ANORM ) / AINVNM

      // Expand AINV into a full matrix and call SSYMM to multiply
      // AINV on the left by A.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 20 J = 1, N
            DO 10 I = 1, J - 1
               AINV( J, I ) = AINV( I, J )
   10       CONTINUE
   20    CONTINUE
      } else {
         DO 40 J = 1, N
            DO 30 I = J + 1, N
               AINV( J, I ) = AINV( I, J )
   30       CONTINUE
   40    CONTINUE
      }
      ssymm('Left', UPLO, N, N, -ONE, A, LDA, AINV, LDAINV, ZERO, WORK, LDWORK );

      // Add the identity matrix to WORK .

      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + ONE
   50 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND ) / EPS ) / REAL( N )

      RETURN

      // End of SPOT03

      }
