      REAL             FUNCTION SQPT01( M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               NORMA
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      SQPT01 = ZERO

      // Test if there is enough workspace

      if ( LWORK.LT.M*N+N ) {
         xerbla('SQPT01', 10 );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK )

      for (J = 1; J <= K; J++) {
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = AF( I, J )
         END DO
         DO I = J + 1, M
            WORK( ( J-1 )*M+I ) = ZERO
         END DO
      END DO
      DO J = K + 1, N
         scopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      END DO

      sormqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare i-th column of QR and jpvt(i)-th column of A

         saxpy(M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      END DO

      SQPT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( REAL( MAX( M, N ) )*SLAMCH( 'Epsilon' ) )       IF( NORMA.NE.ZERO ) SQPT01 = SQPT01 / NORMA

      RETURN

      // End of SQPT01

      }
