      SUBROUTINE SSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ANORM, ULP
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYMM
      // ..
      // .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( N.LE.0 ) RETURN
*
      ULP = SLAMCH( 'Epsilon' )
*
      // Compute product of 1-norms of A and Z.
*
      ANORM = SLANSY( '1', UPLO, N, A, LDA, WORK )* SLANGE( '1', N, M, Z, LDZ, WORK )       IF( ANORM.EQ.ZERO ) ANORM = ONE
*
      IF( ITYPE.EQ.1 ) THEN
*
         // Norm of AZ - BZD
*
         CALL SSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N )
         DO 10 I = 1, M
            CALL SSCAL( N, D( I ), Z( 1, I ), 1 )
   10    CONTINUE
         CALL SSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK, N )
*
         RESULT( 1 ) = ( SLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) / ( N*ULP )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
         // Norm of ABZ - ZD
*
         CALL SSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK, N )
         DO 20 I = 1, M
            CALL SSCAL( N, D( I ), Z( 1, I ), 1 )
   20    CONTINUE
         CALL SSYMM( 'Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z, LDZ )
*
         RESULT( 1 ) = ( SLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
         // Norm of BAZ - ZD
*
         CALL SSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N )
         DO 30 I = 1, M
            CALL SSCAL( N, D( I ), Z( 1, I ), 1 )
   30    CONTINUE
         CALL SSYMM( 'Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z, LDZ )
*
         RESULT( 1 ) = ( SLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP )
      END IF
*
      RETURN
*
      // End of SSGT01
*
      END
