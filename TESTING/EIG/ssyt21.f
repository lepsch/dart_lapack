      SUBROUTINE SSYT21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TEN = 10.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      String             CUPLO;
      int                IINFO, J, JCOL, JR, JROW
      REAL               ANORM, ULP, UNFL, VSAVE, WNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           LSAME, SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLACPY, SLARFY, SLASET, SORM2L, SORM2R, SSYR, SSYR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( ITYPE.EQ.1 ) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         LOWER = .FALSE.
         CUPLO = 'U'
      ELSE
         LOWER = .TRUE.
         CUPLO = 'L'
      END IF
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
*
*     Some Error Checks
*
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         RESULT( 1 ) = TEN / ULP
         RETURN
      END IF
*
*     Do Test 1
*
*     Norm of A:
*
      IF( ITYPE.EQ.3 ) THEN
         ANORM = ONE
      ELSE
         ANORM = MAX( SLANSY( '1', CUPLO, N, A, LDA, WORK ), UNFL )
      END IF
*
*     Compute error matrix:
*
      IF( ITYPE.EQ.1 ) THEN
*
*        ITYPE=1: error = A - U S U**T
*
         CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
         CALL SLACPY( CUPLO, N, N, A, LDA, WORK, N )
*
         DO 10 J = 1, N
            CALL SSYR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK, N )
   10    CONTINUE
*
         IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
            DO 20 J = 1, N - 1
               CALL SSYR2( CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N )
   20       CONTINUE
         END IF
         WNORM = SLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        ITYPE=2: error = V S V**T - A
*
         CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
*
         IF( LOWER ) THEN
            WORK( N**2 ) = D( N )
            DO 40 J = N - 1, 1, -1
               IF( KBAND.EQ.1 ) THEN
                  WORK( ( N+1 )*( J-1 )+2 ) = ( ONE-TAU( J ) )*E( J )
                  DO 30 JR = J + 2, N
                     WORK( ( J-1 )*N+JR ) = -TAU( J )*E( J )*V( JR, J )
   30             CONTINUE
               END IF
*
               VSAVE = V( J+1, J )
               V( J+1, J ) = ONE
               CALL SLARFY( 'L', N-J, V( J+1, J ), 1, TAU( J ), WORK( ( N+1 )*J+1 ), N, WORK( N**2+1 ) )
               V( J+1, J ) = VSAVE
               WORK( ( N+1 )*( J-1 )+1 ) = D( J )
   40       CONTINUE
         ELSE
            WORK( 1 ) = D( 1 )
            DO 60 J = 1, N - 1
               IF( KBAND.EQ.1 ) THEN
                  WORK( ( N+1 )*J ) = ( ONE-TAU( J ) )*E( J )
                  DO 50 JR = 1, J - 1
                     WORK( J*N+JR ) = -TAU( J )*E( J )*V( JR, J+1 )
   50             CONTINUE
               END IF
*
               VSAVE = V( J, J+1 )
               V( J, J+1 ) = ONE
               CALL SLARFY( 'U', J, V( 1, J+1 ), 1, TAU( J ), WORK, N, WORK( N**2+1 ) )
               V( J, J+1 ) = VSAVE
               WORK( ( N+1 )*J+1 ) = D( J+1 )
   60       CONTINUE
         END IF
*
         DO 90 JCOL = 1, N
            IF( LOWER ) THEN
               DO 70 JROW = JCOL, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   70          CONTINUE
            ELSE
               DO 80 JROW = 1, JCOL
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   80          CONTINUE
            END IF
   90    CONTINUE
         WNORM = SLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        ITYPE=3: error = U V**T - I
*
         IF( N.LT.2 ) RETURN
         CALL SLACPY( ' ', N, N, U, LDU, WORK, N )
         IF( LOWER ) THEN
            CALL SORM2R( 'R', 'T', N, N-1, N-1, V( 2, 1 ), LDV, TAU, WORK( N+1 ), N, WORK( N**2+1 ), IINFO )
         ELSE
            CALL SORM2L( 'R', 'T', N, N-1, N-1, V( 1, 2 ), LDV, TAU, WORK, N, WORK( N**2+1 ), IINFO )
         END IF
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = TEN / ULP
            RETURN
         END IF
*
         DO 100 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
  100    CONTINUE
*
         WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
      END IF
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  U U**T - I
*
      IF( ITYPE.EQ.1 ) THEN
         CALL SGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N )
*
         DO 110 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
  110    CONTINUE
*
         RESULT( 2 ) = MIN( SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of SSYT21
*
      END
