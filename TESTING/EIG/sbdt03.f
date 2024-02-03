      SUBROUTINE SBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N;
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                I, J;
      REAL               BNORM, EPS
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SLAMCH
      EXTERNAL           LSAME, ISAMAX, SASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( N.LE.0 ) RETURN
*
*     Compute B - U * S * V' one column at a time.
*
      BNORM = ZERO
      IF( KD.GE.1 ) THEN
*
*        B is bidiagonal.
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           B is upper bidiagonal.
*
            DO 20 J = 1, N
               DO 10 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   10          CONTINUE
               CALL SGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.GT.1 ) THEN
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, SASUM( N, WORK, 1 ) )
   20       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 40 J = 1, N
               DO 30 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   30          CONTINUE
               CALL SGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.LT.N ) THEN
                  WORK( J+1 ) = WORK( J+1 ) + E( J )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, SASUM( N, WORK, 1 ) )
   40       CONTINUE
         END IF
      ELSE
*
*        B is diagonal.
*
         DO 60 J = 1, N
            DO 50 I = 1, N
               WORK( N+I ) = S( I )*VT( I, J )
   50       CONTINUE
            CALL SGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
            WORK( J ) = WORK( J ) + D( J )
            RESID = MAX( RESID, SASUM( N, WORK, 1 ) )
   60    CONTINUE
         J = ISAMAX( N, D, 1 )
         BNORM = ABS( D( J ) )
      END IF
*
*     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
*
      EPS = SLAMCH( 'Precision' )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SBDT03
*
      END
