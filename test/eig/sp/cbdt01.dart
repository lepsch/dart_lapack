      void cbdt01(final int M, final int N, final int KD, final Matrix<double> A_, final int LDA, final Matrix<double> Q_, final int LDQ, final int D, final int E, final Matrix<double> PT_, final int LDPT, final Array<double> _WORK_, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final Q = Q_.dim();
  final PT = PT_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KD, LDA, LDPT, LDQ, M, N;
      double               RESID;
      double               D( * ), E( * ), RWORK( * );
      Complex            A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SCASUM, SLAMCH;
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL

      // Quick return if possible

      if ( M <= 0 || N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute A - Q * B * P**H one column at a time.

      RESID = ZERO;
      if ( KD != 0 ) {

         // B is bidiagonal.

         if ( KD != 0 && M >= N ) {

            // B is upper bidiagonal and M >= N.

            for (J = 1; J <= N; J++) { // 20
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= N - 1; I++) { // 10
                  WORK[M+I] = D( I )*PT( I, J ) + E( I )*PT( I+1, J );
               } // 10
               WORK[M+N] = D( N )*PT( N, J );
               cgemv('No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = max( RESID, SCASUM( M, WORK, 1 ) );
            } // 20
         } else if ( KD < 0 ) {

            // B is upper bidiagonal and M < N.

            for (J = 1; J <= N; J++) { // 40
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= M - 1; I++) { // 30
                  WORK[M+I] = D( I )*PT( I, J ) + E( I )*PT( I+1, J );
               } // 30
               WORK[M+M] = D( M )*PT( M, J );
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = max( RESID, SCASUM( M, WORK, 1 ) );
            } // 40
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 60
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               WORK[M+1] = D( 1 )*PT( 1, J );
               for (I = 2; I <= M; I++) { // 50
                  WORK[M+I] = E( I-1 )*PT( I-1, J ) + D( I )*PT( I, J );
               } // 50
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = max( RESID, SCASUM( M, WORK, 1 ) );
            } // 60
         }
      } else {

         // B is diagonal.

         if ( M >= N ) {
            for (J = 1; J <= N; J++) { // 80
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= N; I++) { // 70
                  WORK[M+I] = D( I )*PT( I, J );
               } // 70
               cgemv('No transpose', M, N, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = max( RESID, SCASUM( M, WORK, 1 ) );
            } // 80
         } else {
            for (J = 1; J <= N; J++) { // 100
               ccopy(M, A( 1, J ), 1, WORK, 1 );
               for (I = 1; I <= M; I++) { // 90
                  WORK[M+I] = D( I )*PT( I, J );
               } // 90
               cgemv('No transpose', M, M, -CMPLX( ONE ), Q, LDQ, WORK( M+1 ), 1, CMPLX( ONE ), WORK, 1 );
               RESID = max( RESID, SCASUM( M, WORK, 1 ) );
            } // 100
         }
      }

      // Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK );
      EPS = SLAMCH( 'Precision' );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS );
         } else {
            if ( ANORM < ONE ) {
               RESID = ( min( RESID, double( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS );
            } else {
               RESID = min( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      }
