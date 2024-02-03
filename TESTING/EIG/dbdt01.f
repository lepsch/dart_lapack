      SUBROUTINE DBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KD, LDA, LDPT, LDQ, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * ), E( * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M.LE.0 .OR. N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute A - Q * B * P**T one column at a time.

      RESID = ZERO
      if ( KD.NE.0 ) {

         // B is bidiagonal.

         if ( KD.NE.0 .AND. M.GE.N ) {

            // B is upper bidiagonal and M >= N.

            DO 20 J = 1, N
               dcopy(M, A( 1, J ), 1, WORK, 1 );
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               dgemv('No transpose', M, N, -ONE, Q, LDQ, WORK( M+1 ), 1, ONE, WORK, 1 );
               RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   20       CONTINUE
         } else if ( KD.LT.0 ) {

            // B is upper bidiagonal and M < N.

            DO 40 J = 1, N
               dcopy(M, A( 1, J ), 1, WORK, 1 );
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK( M+1 ), 1, ONE, WORK, 1 );
               RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   40       CONTINUE
         } else {

            // B is lower bidiagonal.

            DO 60 J = 1, N
               dcopy(M, A( 1, J ), 1, WORK, 1 );
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J )
   50          CONTINUE
               dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK( M+1 ), 1, ONE, WORK, 1 );
               RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   60       CONTINUE
         }
      } else {

         // B is diagonal.

         if ( M.GE.N ) {
            DO 80 J = 1, N
               dcopy(M, A( 1, J ), 1, WORK, 1 );
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               dgemv('No transpose', M, N, -ONE, Q, LDQ, WORK( M+1 ), 1, ONE, WORK, 1 );
               RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   80       CONTINUE
         } else {
            DO 100 J = 1, N
               dcopy(M, A( 1, J ), 1, WORK, 1 );
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK( M+1 ), 1, ONE, WORK, 1 );
               RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
  100       CONTINUE
         }
      }

      // Compute norm(A - Q * B * P**T) / ( n * norm(A) * EPS )

      ANORM = DLANGE( '1', M, N, A, LDA, WORK )
      EPS = DLAMCH( 'Precision' )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         if ( ANORM.GE.RESID ) {
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         } else {
            if ( ANORM.LT.ONE ) {
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS )
            } else {
               RESID = MIN( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            }
         }
      }

      RETURN

      // End of DBDT01

      }
