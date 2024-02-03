      SUBROUTINE DBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J
      DOUBLE PRECISION   BNORM, EPS
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IDAMAX
      DOUBLE PRECISION   DASUM, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DASUM, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
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
               CALL DGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.GT.1 ) THEN
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, DASUM( N, WORK, 1 ) )
   20       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 40 J = 1, N
               DO 30 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   30          CONTINUE
               CALL DGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.LT.N ) THEN
                  WORK( J+1 ) = WORK( J+1 ) + E( J )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, DASUM( N, WORK, 1 ) )
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
            CALL DGEMV( 'No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 )
            WORK( J ) = WORK( J ) + D( J )
            RESID = MAX( RESID, DASUM( N, WORK, 1 ) )
   60    CONTINUE
         J = IDAMAX( N, D, 1 )
         BNORM = ABS( D( J ) )
      END IF
*
*     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
*
      EPS = DLAMCH( 'Precision' )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( DBLE( N )*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, DBLE( N )*BNORM ) / BNORM ) / ( DBLE( N )*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DBDT03
*
      END
