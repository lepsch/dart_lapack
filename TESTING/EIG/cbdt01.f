      SUBROUTINE CBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KD, LDA, LDPT, LDQ, M, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      REAL               CLANGE, SCASUM, SLAMCH
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Compute A - Q * B * P**H one column at a time.

      RESID = ZERO
      IF( KD.NE.0 ) THEN

         // B is bidiagonal.

         IF( KD.NE.0 .AND. M.GE.N ) THEN

            // B is upper bidiagonal and M >= N.

            DO 20 J = 1, N
               CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               CALL CGEMV( 'No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   20       CONTINUE
         ELSE IF( KD.LT.0 ) THEN

            // B is upper bidiagonal and M < N.

            DO 40 J = 1, N
               CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               CALL CGEMV( 'No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   40       CONTINUE
         ELSE

            // B is lower bidiagonal.

            DO 60 J = 1, N
               CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J )
   50          CONTINUE
               CALL CGEMV( 'No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   60       CONTINUE
         END IF
      ELSE

         // B is diagonal.

         IF( M.GE.N ) THEN
            DO 80 J = 1, N
               CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               CALL CGEMV( 'No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   80       CONTINUE
         ELSE
            DO 100 J = 1, N
               CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               CALL CGEMV( 'No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 )
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
  100       CONTINUE
         END IF
      END IF

      // Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      EPS = SLAMCH( 'Precision' )

      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS )
            END IF
         END IF
      END IF

      RETURN

      // End of CBDT01

      END
