      SUBROUTINE ZPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAINV, LDWORK, N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RCOND = ONE
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      AINVNM = ZLANHE( '1', UPLO, N, AINV, LDAINV, RWORK )
      if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      }
      RCOND = ( ONE / ANORM ) / AINVNM

      // Expand AINV into a full matrix and call ZHEMM to multiply
      // AINV on the left by A.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 20 J = 1, N
            DO 10 I = 1, J - 1
               AINV( J, I ) = DCONJG( AINV( I, J ) )
   10       CONTINUE
   20    CONTINUE
      } else {
         DO 40 J = 1, N
            DO 30 I = J + 1, N
               AINV( J, I ) = DCONJG( AINV( I, J ) )
   30       CONTINUE
   40    CONTINUE
      }
      zhemm('Left', UPLO, N, N, -CONE, A, LDA, AINV, LDAINV, CZERO, WORK, LDWORK );

      // Add the identity matrix to WORK .

      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + CONE
   50 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )

      RETURN

      // End of ZPOT03

      }
