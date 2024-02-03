      REAL             FUNCTION SQRT11( M, K, A, LDA, TAU, WORK, LWORK );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Local Arrays ..
      REAL               RDUMMY( 1 );
      // ..
      // .. Executable Statements ..

      SQRT11 = ZERO;

      // Test for sufficient workspace

      if ( LWORK < M*M+M ) {
         xerbla('SQRT11', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0) return;

      slaset('Full', M, M, ZERO, ONE, WORK, M );

      // Form Q

      sorm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      // Form Q'*Q

      sorm2r('Left', 'Transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      for (J = 1; J <= M; J++) {
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE;
      }

      SQRT11 = SLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( REAL( M )*SLAMCH( 'Epsilon' ) );

      return;

      // End of SQRT11

      }
