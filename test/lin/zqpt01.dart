      double zqpt01(M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      Complex         A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             NORMA;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZUNMQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      ZQPT01 = ZERO;

      // Test if there is enough workspace

      if ( LWORK < M*N+N ) {
         xerbla('ZQPT01', 10 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK );

      for (J = 1; J <= K; J++) {
         for (I = 1; I <= min( J, M ); I++) {
            WORK( ( J-1 )*M+I ) = AF( I, J );
         }
         for (I = J + 1; I <= M; I++) {
            WORK( ( J-1 )*M+I ) = ZERO;
         }
      }
      for (J = K + 1; J <= N; J++) {
         zcopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      zunmqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare i-th column of QR and jpvt(i)-th column of A

         zaxpy(M, DCMPLX( -ONE ), A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      ZQPT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( DBLE( max( M, N ) )*DLAMCH( 'Epsilon' ) )       IF( NORMA != ZERO ) ZQPT01 = ZQPT01 / NORMA;

      return;
      }
