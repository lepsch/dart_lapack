      double           FUNCTION ZRZT01( M, N, A, AF, LDA, TAU, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             NORMA;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZLASET, ZUNMRZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      ZRZT01 = ZERO

      if ( LWORK.LT.M*N+M ) {
         xerbla('ZRZT01', 8 );
         RETURN
      }

      // Quick return if possible

      if (M.LE.0 .OR. N.LE.0) RETURN;

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )

      // Copy upper triangle R

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK, M );
      for (J = 1; J <= M; J++) { // 20
         for (I = 1; I <= J; I++) { // 10
            WORK( ( J-1 )*M+I ) = AF( I, J )
         } // 10
      } // 20

      // R = R * P(1) * ... *P(m)

      zunmrz('Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      // R = R - A

      for (I = 1; I <= N; I++) { // 30
         zaxpy(M, DCMPLX( -ONE ), A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 );
      } // 30

      ZRZT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK )

      ZRZT01 = ZRZT01 / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      if (NORMA.NE.ZERO) ZRZT01 = ZRZT01 / NORMA;

      RETURN

      // End of ZRZT01

      }
