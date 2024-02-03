      SUBROUTINE CGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      REAL               RESULT
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      REAL               ANORM, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      RESULT = ZERO
      IF( N.LE.0 ) RETURN
*
*     Constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
*
*     Some Error Checks
*
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         RESULT = TEN / ULP
         RETURN
      END IF
*
      IF( ITYPE.LE.2 ) THEN
*
*        Tests scaled by the norm(A)
*
         ANORM = MAX( CLANGE( '1', N, N, A, LDA, RWORK ), UNFL )
*
         IF( ITYPE.EQ.1 ) THEN
*
*           ITYPE=1: Compute W = A - U B V**H
*
            CALL CLACPY( ' ', N, N, A, LDA, WORK, N )
            CALL CGEMM( 'N', 'N', N, N, N, CONE, U, LDU, B, LDB, CZERO, WORK( N**2+1 ), N )
*
            CALL CGEMM( 'N', 'C', N, N, N, -CONE, WORK( N**2+1 ), N, V, LDV, CONE, WORK, N )
*
         ELSE
*
*           ITYPE=2: Compute W = A - B
*
            CALL CLACPY( ' ', N, N, B, LDB, WORK, N )
*
            DO 20 JCOL = 1, N
               DO 10 JROW = 1, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         END IF
*
*        Compute norm(W)/ ( ulp*norm(A) )
*
         WNORM = CLANGE( '1', N, N, WORK, N, RWORK )
*
         IF( ANORM.GT.WNORM ) THEN
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            ELSE
               RESULT = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
            END IF
         END IF
*
      ELSE
*
*        Tests not scaled by norm(A)
*
*        ITYPE=3: Compute  U U**H - I
*
         CALL CGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N )
*
         DO 30 JDIAG = 1, N
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - CONE
   30    CONTINUE
*
         RESULT = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of CGET51
*
      END
