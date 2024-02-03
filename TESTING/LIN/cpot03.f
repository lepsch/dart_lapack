      SUBROUTINE CPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAINV, LDWORK, N;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               AINVNM, ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL LSAME, CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
      AINVNM = CLANHE( '1', UPLO, N, AINV, LDAINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE/ANORM ) / AINVNM

      // Expand AINV into a full matrix and call CHEMM to multiply
      // AINV on the left by A.

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J - 1
               AINV( J, I ) = CONJG( AINV( I, J ) )
   10       CONTINUE
   20    CONTINUE
      } else {
         DO 40 J = 1, N
            DO 30 I = J + 1, N
               AINV( J, I ) = CONJG( AINV( I, J ) )
   30       CONTINUE
   40    CONTINUE
      END IF
      CALL CHEMM( 'Left', UPLO, N, N, -CONE, A, LDA, AINV, LDAINV, CZERO, WORK, LDWORK )

      // Add the identity matrix to WORK .

      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + CONE
   50 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N )

      RETURN

      // End of CPOT03

      }
