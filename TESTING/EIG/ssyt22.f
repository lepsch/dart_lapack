      SUBROUTINE SSYT22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                J, JJ, JJ1, JJ2, NN, NNP1
      REAL               ANORM, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANSY
      EXTERNAL           SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SSYMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 .OR. M.LE.0 ) RETURN
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )
*
*     Do Test 1
*
*     Norm of A:
*
      ANORM = MAX( SLANSY( '1', UPLO, N, A, LDA, WORK ), UNFL )
*
*     Compute error matrix:
*
*     ITYPE=1: error = U**T A U - S
*
      CALL SSYMM( 'L', UPLO, N, M, ONE, A, LDA, U, LDU, ZERO, WORK, N )
      NN = N*N
      NNP1 = NN + 1
      CALL SGEMM( 'T', 'N', M, M, N, ONE, U, LDU, WORK, N, ZERO, WORK( NNP1 ), N )
      DO 10 J = 1, M
         JJ = NN + ( J-1 )*N + J
         WORK( JJ ) = WORK( JJ ) - D( J )
   10 CONTINUE
      IF( KBAND.EQ.1 .AND. N.GT.1 ) THEN
         DO 20 J = 2, M
            JJ1 = NN + ( J-1 )*N + J - 1
            JJ2 = NN + ( J-2 )*N + J
            WORK( JJ1 ) = WORK( JJ1 ) - E( J-1 )
            WORK( JJ2 ) = WORK( JJ2 ) - E( J-1 )
   20    CONTINUE
      END IF
      WNORM = SLANSY( '1', UPLO, M, WORK( NNP1 ), N, WORK( 1 ) )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( M ) ) / ( M*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  U**T U - I
*
      IF( ITYPE.EQ.1 ) CALL SORT01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RESULT( 2 ) )
*
      RETURN
*
*     End of SSYT22
*
      END
