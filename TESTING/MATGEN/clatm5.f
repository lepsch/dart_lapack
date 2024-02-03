      SUBROUTINE CLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA, QBLCKB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N, PRTYPE, QBLCKA, QBLCKB;
      REAL               ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * ), L( LDL, * ), R( LDR, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, TWO, ZERO, HALF, TWENTY
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ), TWO = ( 2.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ), TWENTY = ( 2.0E+1, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, J, K;
      COMPLEX            IMEPS, REEPS
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MOD, SIN
*     ..
*     .. External Subroutines ..
      // EXTERNAL CGEMM
*     ..
*     .. Executable Statements ..
*
      IF( PRTYPE.EQ.1 ) THEN
         DO 20 I = 1, M
            DO 10 J = 1, M
               IF( I.EQ.J ) THEN
                  A( I, J ) = ONE
                  D( I, J ) = ONE
               ELSE IF( I.EQ.J-1 ) THEN
                  A( I, J ) = -ONE
                  D( I, J ) = ZERO
               ELSE
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               END IF
   10       CONTINUE
   20    CONTINUE
*
         DO 40 I = 1, N
            DO 30 J = 1, N
               IF( I.EQ.J ) THEN
                  B( I, J ) = ONE - ALPHA
                  E( I, J ) = ONE
               ELSE IF( I.EQ.J-1 ) THEN
                  B( I, J ) = ONE
                  E( I, J ) = ZERO
               ELSE
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               END IF
   30       CONTINUE
   40    CONTINUE
*
         DO 60 I = 1, M
            DO 50 J = 1, N
               R( I, J ) = ( HALF-SIN( CMPLX( I / J ) ) )*TWENTY
               L( I, J ) = R( I, J )
   50       CONTINUE
   60    CONTINUE
*
      ELSE IF( PRTYPE.EQ.2 .OR. PRTYPE.EQ.3 ) THEN
         DO 80 I = 1, M
            DO 70 J = 1, M
               IF( I.LE.J ) THEN
                  A( I, J ) = ( HALF-SIN( CMPLX( I ) ) )*TWO
                  D( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
               ELSE
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               END IF
   70       CONTINUE
   80    CONTINUE
*
         DO 100 I = 1, N
            DO 90 J = 1, N
               IF( I.LE.J ) THEN
                  B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
                  E( I, J ) = ( HALF-SIN( CMPLX( J ) ) )*TWO
               ELSE
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               END IF
   90       CONTINUE
  100    CONTINUE
*
         DO 120 I = 1, M
            DO 110 J = 1, N
               R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
  110       CONTINUE
  120    CONTINUE
*
         IF( PRTYPE.EQ.3 ) THEN
            IF( QBLCKA.LE.1 ) QBLCKA = 2
            DO 130 K = 1, M - 1, QBLCKA
               A( K+1, K+1 ) = A( K, K )
               A( K+1, K ) = -SIN( A( K, K+1 ) )
  130       CONTINUE
*
            IF( QBLCKB.LE.1 ) QBLCKB = 2
            DO 140 K = 1, N - 1, QBLCKB
               B( K+1, K+1 ) = B( K, K )
               B( K+1, K ) = -SIN( B( K, K+1 ) )
  140       CONTINUE
         END IF
*
      ELSE IF( PRTYPE.EQ.4 ) THEN
         DO 160 I = 1, M
            DO 150 J = 1, M
               A( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
               D( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
  150       CONTINUE
  160    CONTINUE
*
         DO 180 I = 1, N
            DO 170 J = 1, N
               B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
               E( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
  170       CONTINUE
  180    CONTINUE
*
         DO 200 I = 1, M
            DO 190 J = 1, N
               R( I, J ) = ( HALF-SIN( CMPLX( J / I ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
  190       CONTINUE
  200    CONTINUE
*
      ELSE IF( PRTYPE.GE.5 ) THEN
         REEPS = HALF*TWO*TWENTY / ALPHA
         IMEPS = ( HALF-TWO ) / ALPHA
         DO 220 I = 1, M
            DO 210 J = 1, N
               R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*ALPHA / TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*ALPHA / TWENTY
  210       CONTINUE
  220    CONTINUE
*
         DO 230 I = 1, M
            D( I, I ) = ONE
  230    CONTINUE
*
         DO 240 I = 1, M
            IF( I.LE.4 ) THEN
               A( I, I ) = ONE
               IF( I.GT.2 ) A( I, I ) = ONE + REEPS
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.M ) THEN
                  A( I, I+1 ) = IMEPS
               ELSE IF( I.GT.1 ) THEN
                  A( I, I-1 ) = -IMEPS
               END IF
            ELSE IF( I.LE.8 ) THEN
               IF( I.LE.6 ) THEN
                  A( I, I ) = REEPS
               ELSE
                  A( I, I ) = -REEPS
               END IF
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.M ) THEN
                  A( I, I+1 ) = ONE
               ELSE IF( I.GT.1 ) THEN
                  A( I, I-1 ) = -ONE
               END IF
            ELSE
               A( I, I ) = ONE
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.M ) THEN
                  A( I, I+1 ) = IMEPS*2
               ELSE IF( I.GT.1 ) THEN
                  A( I, I-1 ) = -IMEPS*2
               END IF
            END IF
  240    CONTINUE
*
         DO 250 I = 1, N
            E( I, I ) = ONE
            IF( I.LE.4 ) THEN
               B( I, I ) = -ONE
               IF( I.GT.2 ) B( I, I ) = ONE - REEPS
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.N ) THEN
                  B( I, I+1 ) = IMEPS
               ELSE IF( I.GT.1 ) THEN
                  B( I, I-1 ) = -IMEPS
               END IF
            ELSE IF( I.LE.8 ) THEN
               IF( I.LE.6 ) THEN
                  B( I, I ) = REEPS
               ELSE
                  B( I, I ) = -REEPS
               END IF
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.N ) THEN
                  B( I, I+1 ) = ONE + IMEPS
               ELSE IF( I.GT.1 ) THEN
                  B( I, I-1 ) = -ONE - IMEPS
               END IF
            ELSE
               B( I, I ) = ONE - REEPS
               IF( MOD( I, 2 ).NE.0 .AND. I.LT.N ) THEN
                  B( I, I+1 ) = IMEPS*2
               ELSE IF( I.GT.1 ) THEN
                  B( I, I-1 ) = -IMEPS*2
               END IF
            END IF
  250    CONTINUE
      END IF
*
*     Compute rhs (C, F)
*
      CALL CGEMM( 'N', 'N', M, N, M, ONE, A, LDA, R, LDR, ZERO, C, LDC )
      CALL CGEMM( 'N', 'N', M, N, N, -ONE, L, LDL, B, LDB, ONE, C, LDC )
      CALL CGEMM( 'N', 'N', M, N, M, ONE, D, LDD, R, LDR, ZERO, F, LDF )
      CALL CGEMM( 'N', 'N', M, N, N, -ONE, L, LDL, E, LDE, ONE, F, LDF )
*
*     End of CLATM5
*
      END
