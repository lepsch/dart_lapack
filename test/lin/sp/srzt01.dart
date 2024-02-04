      double srzt01(M, N, A, AF, LDA, TAU, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double               NORMA;
      // ..
      // .. Local Arrays ..
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SLASET, SORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..

      SRZT01 = ZERO;

      if ( LWORK < M*N+M ) {
         xerbla('SRZT01', 8 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK );

      // Copy upper triangle R

      slaset('Full', M, N, ZERO, ZERO, WORK, M );
      for (J = 1; J <= M; J++) { // 20
         for (I = 1; I <= J; I++) { // 10
            WORK[( J-1 )*M+I] = AF( I, J );
         } // 10
      } // 20

      // R = R * P(1) * ... *P(m)

      sormrz('Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      // R = R - A

      for (I = 1; I <= N; I++) { // 30
         saxpy(M, -ONE, A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 );
      } // 30

      SRZT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK );

      SRZT01 = SRZT01 / ( SLAMCH( 'Epsilon' )*REAL( max( M, N ) ) );
      if (NORMA != ZERO) SRZT01 = SRZT01 / NORMA;

      return;
      }
