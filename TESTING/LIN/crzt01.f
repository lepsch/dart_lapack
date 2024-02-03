      REAL             FUNCTION CRZT01( M, N, A, AF, LDA, TAU, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               NORMA;
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 );
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CLASET, CUNMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Executable Statements ..

      CRZT01 = ZERO;

      if ( LWORK < M*N+M ) {
         xerbla('CRZT01', 8 );
         RETURN;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) RETURN;

      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK );

      // Copy upper triangle R

      claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), WORK, M );
      for (J = 1; J <= M; J++) { // 20
         for (I = 1; I <= J; I++) { // 10
            WORK( ( J-1 )*M+I ) = AF( I, J );
         } // 10
      } // 20

      // R = R * P(1) * ... *P(m)

      cunmrz('Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      // R = R - A

      for (I = 1; I <= N; I++) { // 30
         caxpy(M, CMPLX( -ONE ), A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 );
      } // 30

      CRZT01 = CLANGE( 'One-norm', M, N, WORK, M, RWORK );

      CRZT01 = CRZT01 / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) );
      if (NORMA != ZERO) CRZT01 = CRZT01 / NORMA;

      RETURN;

      // End of CRZT01

      }
