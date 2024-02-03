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
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
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

      if ( M.LE.0 .OR. N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute A - Q * B * P**H one column at a time.

      RESID = ZERO
      if ( KD.NE.0 ) {

         // B is bidiagonal.

         if ( KD.NE.0 .AND. M.GE.N ) {

            // B is upper bidiagonal and M >= N.

            DO 20 J = 1, N
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               cgemv('No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   20       CONTINUE
         } else if ( KD.LT.0 ) {

            // B is upper bidiagonal and M < N.

            DO 40 J = 1, N
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   40       CONTINUE
         } else {

            // B is lower bidiagonal.

            DO 60 J = 1, N
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J )
   50          CONTINUE
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   60       CONTINUE
         }
      } else {

         // B is diagonal.

         if ( M.GE.N ) {
            DO 80 J = 1, N
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               cgemv('No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   80       CONTINUE
         } else {
            DO 100 J = 1, N
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
  100       CONTINUE
         }
      }

      // Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      EPS = SLAMCH( 'Precision' )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         if ( ANORM.GE.RESID ) {
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
         } else {
            if ( ANORM.LT.ONE ) {
               RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS )
            } else {
               RESID = MIN( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS )
            }
         }
      }

      RETURN

      // End of CBDT01

      }
