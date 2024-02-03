      SUBROUTINE ZSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, ULP;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZHEMM
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      IF( N.LE.0 ) RETURN

      ULP = DLAMCH( 'Epsilon' )

      // Compute product of 1-norms of A and Z.

      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )* ZLANGE( '1', N, M, Z, LDZ, RWORK )       IF( ANORM.EQ.ZERO ) ANORM = ONE

      IF( ITYPE.EQ.1 ) THEN

         // Norm of AZ - BZD

         CALL ZHEMM( 'Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N )
         DO 10 I = 1, M
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   10    CONTINUE
         CALL ZHEMM( 'Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, -CONE, WORK, N )

         RESULT( 1 ) = ( ZLANGE( '1', N, M, WORK, N, RWORK ) / ANORM ) / ( N*ULP )

      ELSE IF( ITYPE.EQ.2 ) THEN

         // Norm of ABZ - ZD

         CALL ZHEMM( 'Left', UPLO, N, M, CONE, B, LDB, Z, LDZ, CZERO, WORK, N )
         DO 20 I = 1, M
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   20    CONTINUE
         CALL ZHEMM( 'Left', UPLO, N, M, CONE, A, LDA, WORK, N, -CONE, Z, LDZ )

         RESULT( 1 ) = ( ZLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )

      ELSE IF( ITYPE.EQ.3 ) THEN

         // Norm of BAZ - ZD

         CALL ZHEMM( 'Left', UPLO, N, M, CONE, A, LDA, Z, LDZ, CZERO, WORK, N )
         DO 30 I = 1, M
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   30    CONTINUE
         CALL ZHEMM( 'Left', UPLO, N, M, CONE, B, LDB, WORK, N, -CONE, Z, LDZ )

         RESULT( 1 ) = ( ZLANGE( '1', N, M, Z, LDZ, RWORK ) / ANORM ) / ( N*ULP )
      END IF

      RETURN

      // End of ZDGT01

      END
