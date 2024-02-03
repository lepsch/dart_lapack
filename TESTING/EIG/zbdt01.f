      SUBROUTINE ZBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                KD, LDA, LDPT, LDQ, M, N;
      double             RESID;
*     ..
*     .. Array Arguments ..
      double             D( * ), E( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
*     ..
*     .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      EXTERNAL           DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Compute A - Q * B * P**H one column at a time.
*
      RESID = ZERO
      IF( KD.NE.0 ) THEN
*
*        B is bidiagonal.
*
         IF( KD.NE.0 .AND. M.GE.N ) THEN
*
*           B is upper bidiagonal and M >= N.
*
            DO 20 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               CALL ZGEMV( 'No transpose', M, N, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   20       CONTINUE
         ELSE IF( KD.LT.0 ) THEN
*
*           B is upper bidiagonal and M < N.
*
            DO 40 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   40       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 60 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J )
   50          CONTINUE
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   60       CONTINUE
         END IF
      ELSE
*
*        B is diagonal.
*
         IF( M.GE.N ) THEN
            DO 80 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               CALL ZGEMV( 'No transpose', M, N, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   80       CONTINUE
         ELSE
            DO 100 J = 1, N
               CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
  100       CONTINUE
         END IF
      END IF
*
*     Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      EPS = DLAMCH( 'Precision' )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZBDT01
*
      END
