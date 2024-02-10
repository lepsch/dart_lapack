      double drzt01(M, N, A, AF, LDA, TAU, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J;
      double             NORMA;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLASET, DORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX

      DRZT01 = ZERO;

      if ( LWORK < M*N+M ) {
         xerbla('DRZT01', 8 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = dlange( 'One-norm', M, N, A, LDA, RWORK );

      // Copy upper triangle R

      dlaset('Full', M, N, ZERO, ZERO, WORK, M );
      for (J = 1; J <= M; J++) { // 20
         for (I = 1; I <= J; I++) { // 10
            WORK[( J-1 )*M+I] = AF( I, J );
         } // 10
      } // 20

      // R = R * P(1) * ... *P(m)

      dormrz('Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      // R = R - A

      for (I = 1; I <= N; I++) { // 30
         daxpy(M, -ONE, A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 );
      } // 30

      DRZT01 = dlange( 'One-norm', M, N, WORK, M, RWORK );

      DRZT01 = DRZT01 / ( dlamch( 'Epsilon' )*(max( M, N )).toDouble() );
      if (NORMA != ZERO) DRZT01 = DRZT01 / NORMA;

      }
