      SUBROUTINE SSPT21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP, TAU, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDU, N
*     ..
*     .. Array Arguments ..
      REAL               AP( * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), VP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TEN = 10.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 1.0E+0 / 2.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JP, JP1, JR, LAP
      REAL               ANORM, TEMP, ULP, UNFL, VSAVE, WNORM
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SDOT, SLAMCH, SLANGE, SLANSP
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANGE, SLANSP
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMM, SLACPY, SLASET, SOPMTR, SSPMV, SSPR, SSPR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     1)      Constants
*
      RESULT( 1 ) = ZERO
      IF( ITYPE.EQ.1 ) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN
*
      LAP = ( N*( N+1 ) ) / 2
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
         ANORM = MAX( SLANSP( '1', CUPLO, N, AP, WORK ), UNFL )
      END IF
*
*     Compute error matrix:
*
      IF( ITYPE.EQ.1 ) THEN
*
*        ITYPE=1: error = A - U S U**T
*
         CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
         CALL SCOPY( LAP, AP, 1, WORK, 1 )
*
         DO 10 J = 1, N
            CALL SSPR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK )
   10    CONTINUE
*
         IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
            DO 20 J = 1, N - 1
               CALL SSPR2( CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK )
   20       CONTINUE
         END IF
         WNORM = SLANSP( '1', CUPLO, N, WORK, WORK( N**2+1 ) )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        ITYPE=2: error = V S V**T - A
*
         CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
*
         IF( LOWER ) THEN
            WORK( LAP ) = D( N )
            DO 40 J = N - 1, 1, -1
               JP = ( ( 2*N-J )*( J-1 ) ) / 2
               JP1 = JP + N - J
               IF( KBAND.EQ.1 ) THEN
                  WORK( JP+J+1 ) = ( ONE-TAU( J ) )*E( J )
                  DO 30 JR = J + 2, N
                     WORK( JP+JR ) = -TAU( J )*E( J )*VP( JP+JR )
   30             CONTINUE
               END IF
*
               IF( TAU( J ).NE.ZERO ) THEN
                  VSAVE = VP( JP+J+1 )
                  VP( JP+J+1 ) = ONE
                  CALL SSPMV( 'L', N-J, ONE, WORK( JP1+J+1 ), VP( JP+J+1 ), 1, ZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*SDOT( N-J, WORK( LAP+1 ), 1, VP( JP+J+1 ), 1 )                   CALL SAXPY( N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), 1 )                   CALL SSPR2( 'L', N-J, -TAU( J ), VP( JP+J+1 ), 1, WORK( LAP+1 ), 1, WORK( JP1+J+1 ) )
                  VP( JP+J+1 ) = VSAVE
               END IF
               WORK( JP+J ) = D( J )
   40       CONTINUE
         ELSE
            WORK( 1 ) = D( 1 )
            DO 60 J = 1, N - 1
               JP = ( J*( J-1 ) ) / 2
               JP1 = JP + J
               IF( KBAND.EQ.1 ) THEN
                  WORK( JP1+J ) = ( ONE-TAU( J ) )*E( J )
                  DO 50 JR = 1, J - 1
                     WORK( JP1+JR ) = -TAU( J )*E( J )*VP( JP1+JR )
   50             CONTINUE
               END IF
*
               IF( TAU( J ).NE.ZERO ) THEN
                  VSAVE = VP( JP1+J )
                  VP( JP1+J ) = ONE
                  CALL SSPMV( 'U', J, ONE, WORK, VP( JP1+1 ), 1, ZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*SDOT( J, WORK( LAP+1 ), 1, VP( JP1+1 ), 1 )                   CALL SAXPY( J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), 1 )                   CALL SSPR2( 'U', J, -TAU( J ), VP( JP1+1 ), 1, WORK( LAP+1 ), 1, WORK )
                  VP( JP1+J ) = VSAVE
               END IF
               WORK( JP1+J+1 ) = D( J+1 )
   60       CONTINUE
         END IF
*
         DO 70 J = 1, LAP
            WORK( J ) = WORK( J ) - AP( J )
   70    CONTINUE
         WNORM = SLANSP( '1', CUPLO, N, WORK, WORK( LAP+1 ) )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        ITYPE=3: error = U V**T - I
*
         IF( N.LT.2 ) RETURN
         CALL SLACPY( ' ', N, N, U, LDU, WORK, N )
         CALL SOPMTR( 'R', CUPLO, 'T', N, N, VP, TAU, WORK, N, WORK( N**2+1 ), IINFO )
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = TEN / ULP
            RETURN
         END IF
*
         DO 80 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   80    CONTINUE
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
         DO 90 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   90    CONTINUE
*
         RESULT( 2 ) = MIN( SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of SSPT21
*
      END
