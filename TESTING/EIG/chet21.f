      SUBROUTINE CHET21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JCOL, JR, JROW;
      REAL               ANORM, ULP, UNFL, WNORM
      COMPLEX            VSAVE
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL LSAME, CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHER, CHER2, CLACPY, CLARFY, CLASET, CUNM2L, CUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..
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
      // Some Error Checks
*
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         RESULT( 1 ) = TEN / ULP
         RETURN
      END IF
*
      // Do Test 1
*
      // Norm of A:
*
      IF( ITYPE.EQ.3 ) THEN
         ANORM = ONE
      ELSE
         ANORM = MAX( CLANHE( '1', CUPLO, N, A, LDA, RWORK ), UNFL )
      END IF
*
      // Compute error matrix:
*
      IF( ITYPE.EQ.1 ) THEN
*
         // ITYPE=1: error = A - U S U**H
*
         CALL CLASET( 'Full', N, N, CZERO, CZERO, WORK, N )
         CALL CLACPY( CUPLO, N, N, A, LDA, WORK, N )
*
         DO 10 J = 1, N
            CALL CHER( CUPLO, N, -D( J ), U( 1, J ), 1, WORK, N )
   10    CONTINUE
*
         IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
            DO 20 J = 1, N - 1
               CALL CHER2( CUPLO, N, -CMPLX( E( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N )
   20       CONTINUE
         END IF
         WNORM = CLANHE( '1', CUPLO, N, WORK, N, RWORK )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
         // ITYPE=2: error = V S V**H - A
*
         CALL CLASET( 'Full', N, N, CZERO, CZERO, WORK, N )
*
         IF( LOWER ) THEN
            WORK( N**2 ) = D( N )
            DO 40 J = N - 1, 1, -1
               IF( KBAND.EQ.1 ) THEN
                  WORK( ( N+1 )*( J-1 )+2 ) = ( CONE-TAU( J ) )*E( J )
                  DO 30 JR = J + 2, N
                     WORK( ( J-1 )*N+JR ) = -TAU( J )*E( J )*V( JR, J )
   30             CONTINUE
               END IF
*
               VSAVE = V( J+1, J )
               V( J+1, J ) = ONE
               CALL CLARFY( 'L', N-J, V( J+1, J ), 1, TAU( J ), WORK( ( N+1 )*J+1 ), N, WORK( N**2+1 ) )
               V( J+1, J ) = VSAVE
               WORK( ( N+1 )*( J-1 )+1 ) = D( J )
   40       CONTINUE
         ELSE
            WORK( 1 ) = D( 1 )
            DO 60 J = 1, N - 1
               IF( KBAND.EQ.1 ) THEN
                  WORK( ( N+1 )*J ) = ( CONE-TAU( J ) )*E( J )
                  DO 50 JR = 1, J - 1
                     WORK( J*N+JR ) = -TAU( J )*E( J )*V( JR, J+1 )
   50             CONTINUE
               END IF
*
               VSAVE = V( J, J+1 )
               V( J, J+1 ) = ONE
               CALL CLARFY( 'U', J, V( 1, J+1 ), 1, TAU( J ), WORK, N, WORK( N**2+1 ) )
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
         WNORM = CLANHE( '1', CUPLO, N, WORK, N, RWORK )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
         // ITYPE=3: error = U V**H - I
*
         IF( N.LT.2 ) RETURN
         CALL CLACPY( ' ', N, N, U, LDU, WORK, N )
         IF( LOWER ) THEN
            CALL CUNM2R( 'R', 'C', N, N-1, N-1, V( 2, 1 ), LDV, TAU, WORK( N+1 ), N, WORK( N**2+1 ), IINFO )
         ELSE
            CALL CUNM2L( 'R', 'C', N, N-1, N-1, V( 1, 2 ), LDV, TAU, WORK, N, WORK( N**2+1 ), IINFO )
         END IF
         IF( IINFO.NE.0 ) THEN
            RESULT( 1 ) = TEN / ULP
            RETURN
         END IF
*
         DO 100 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
  100    CONTINUE
*
         WNORM = CLANGE( '1', N, N, WORK, N, RWORK )
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
      // Do Test 2
*
      // Compute  U U**H - I
*
      IF( ITYPE.EQ.1 ) THEN
         CALL CGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N )
*
         DO 110 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
  110    CONTINUE
*
         RESULT( 2 ) = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP )
      END IF
*
      RETURN
*
      // End of CHET21
*
      END
