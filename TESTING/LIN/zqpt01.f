      double           FUNCTION ZQPT01( M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
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
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZUNMQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      ZQPT01 = ZERO

      // Test if there is enough workspace

      if ( LWORK.LT.M*N+N ) {
         xerbla('ZQPT01', 10 );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )

      for (J = 1; J <= K; J++) {
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = AF( I, J )
         END DO
         DO I = J + 1, M
            WORK( ( J-1 )*M+I ) = ZERO
         END DO
      END DO
      DO J = K + 1, N
         zcopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      END DO

      zunmqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare i-th column of QR and jpvt(i)-th column of A

         zaxpy(M, DCMPLX( -ONE ), A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      END DO

      ZQPT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( DBLE( MAX( M, N ) )*DLAMCH( 'Epsilon' ) )       IF( NORMA.NE.ZERO ) ZQPT01 = ZQPT01 / NORMA

      RETURN

      // End of ZQPT01

      }
