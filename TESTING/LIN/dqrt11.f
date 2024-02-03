      double           FUNCTION DQRT11( M, K, A, LDA, TAU, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Local Arrays ..
      double             RDUMMY( 1 );
      // ..
      // .. Executable Statements ..

      DQRT11 = ZERO;

      // Test for sufficient workspace

      if ( LWORK < M*M+M ) {
         xerbla('DQRT11', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0) RETURN;

      dlaset('Full', M, M, ZERO, ONE, WORK, M );

      // Form Q

      dorm2r('Left', 'No transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      // Form Q'*Q

      dorm2r('Left', 'Transpose', M, M, K, A, LDA, TAU, WORK, M, WORK( M*M+1 ), INFO );

      for (J = 1; J <= M; J++) {
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE;
      }

      DQRT11 = DLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) / ( DBLE( M )*DLAMCH( 'Epsilon' ) );

      return;

      // End of DQRT11

      }
