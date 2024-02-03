      SUBROUTINE ZBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KD, LDA, LDPT, LDQ, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESID = ZERO
         RETURN
      }

      // Compute A - Q * B * P**H one column at a time.

      RESID = ZERO
      if ( KD != 0 ) {

         // B is bidiagonal.

         if ( KD != 0 && M >= N ) {

            // B is upper bidiagonal and M >= N.

            for (J = 1; J <= N; J++) { // 20
               zcopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= N - 1; I++) { // 10
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
               } // 10
               WORK( M+N ) = D( N )*PT( N, J )
               zgemv('No transpose', M, N, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
            } // 20
         } else if ( KD < 0 ) {

            // B is upper bidiagonal and M < N.

            for (J = 1; J <= N; J++) { // 40
               zcopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= M - 1; I++) { // 30
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
               } // 30
               WORK( M+M ) = D( M )*PT( M, J )
               zgemv('No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
            } // 40
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 60
               zcopy(M, A( 1, J ), 1, WORK, 1 );
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               for (I = 2; I <= M; I++) { // 50
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J )
               } // 50
               zgemv('No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
            } // 60
         }
      } else {

         // B is diagonal.

         if ( M >= N ) {
            for (J = 1; J <= N; J++) { // 80
               zcopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= N; I++) { // 70
                  WORK( M+I ) = D( I )*PT( I, J )
               } // 70
               zgemv('No transpose', M, N, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
            } // 80
         } else {
            for (J = 1; J <= N; J++) { // 100
               zcopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= M; I++) { // 90
                  WORK( M+I ) = D( I )*PT( I, J )
               } // 90
               zgemv('No transpose', M, M, -DCMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, DCMPLX( ONE ), WORK, 1 );
               RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
            } // 100
         }
      }

      // Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )

      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      EPS = DLAMCH( 'Precision' )

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         } else {
            if ( ANORM < ONE ) {
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS )
            } else {
               RESID = MIN( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            }
         }
      }

      RETURN

      // End of ZBDT01

      }
