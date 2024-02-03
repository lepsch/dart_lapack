      SUBROUTINE ZHET22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
*     ..
*     .. Array Arguments ..
      double             D( * ), E( * ), RESULT( 2 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      double             ANORM, ULP, UNFL, WNORM;
*     ..
*     .. External Functions ..
      double             DLAMCH, ZLANHE;
      // EXTERNAL DLAMCH, ZLANHE
*     ..
*     .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHEMM, ZUNT01
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 .OR. M.LE.0 ) RETURN
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
*
*     Do Test 1
*
*     Norm of A:
*
      ANORM = MAX( ZLANHE( '1', UPLO, N, A, LDA, RWORK ), UNFL )
*
*     Compute error matrix:
*
*     ITYPE=1: error = U**H A U - S
*
      CALL ZHEMM( 'L', UPLO, N, M, CONE, A, LDA, U, LDU, CZERO, WORK, N )
      NN = N*N
      NNP1 = NN + 1
      CALL ZGEMM( 'C', 'N', M, M, N, CONE, U, LDU, WORK, N, CZERO, WORK( NNP1 ), N )
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
      WNORM = ZLANHE( '1', UPLO, M, WORK( NNP1 ), N, RWORK )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  U**H U - I
*
      IF( ITYPE.EQ.1 ) CALL ZUNT01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RWORK, RESULT( 2 ) )
*
      RETURN
*
*     End of ZHET22
*
      END
