      SUBROUTINE DSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, LDA, LDB, LDZ, M, N
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      int                I
      double             ANORM, ULP;
*     ..
*     .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSYMM
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( N.LE.0 ) RETURN
*
      ULP = DLAMCH( 'Epsilon' )
*
*     Compute product of 1-norms of A and Z.
*
      ANORM = DLANSY( '1', UPLO, N, A, LDA, WORK )* DLANGE( '1', N, M, Z, LDZ, WORK )       IF( ANORM.EQ.ZERO ) ANORM = ONE
*
      IF( ITYPE.EQ.1 ) THEN
*
*        Norm of AZ - BZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N )
         DO 10 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   10    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK, N )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) / ( N*ULP )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Norm of ABZ - ZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK, N )
         DO 20 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   20    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z, LDZ )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Norm of BAZ - ZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK, N )
         DO 30 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   30    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z, LDZ )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of DSGT01
*
      END
