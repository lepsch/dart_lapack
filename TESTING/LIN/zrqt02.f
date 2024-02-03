      SUBROUTINE ZRQT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RESULT( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         ROGUE
      PARAMETER          ( ROGUE = ( -1.0D+10, -1.0D+10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGRQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the last k rows of the factorization to the array Q
*
      CALL ZLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      IF( K.LT.N ) CALL ZLACPY( 'Full', K, N-K, AF( M-K+1, 1 ), LDA, Q( M-K+1, 1 ), LDA )       IF( K.GT.1 ) CALL ZLACPY( 'Lower', K-1, K-1, AF( M-K+2, N-K+1 ), LDA, Q( M-K+2, N-K+1 ), LDA )
*
*     Generate the last n rows of the matrix Q
*
      SRNAMT = 'ZUNGRQ'
      CALL ZUNGRQ( M, N, K, Q, LDA, TAU( M-K+1 ), WORK, LWORK, INFO )
*
*     Copy R(m-k+1:m,n-m+1:n)
*
      CALL ZLASET( 'Full', K, M, DCMPLX( ZERO ), DCMPLX( ZERO ), R( M-K+1, N-M+1 ), LDA )       CALL ZLACPY( 'Upper', K, K, AF( M-K+1, N-K+1 ), LDA, R( M-K+1, N-K+1 ), LDA )
*
*     Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'
*
      CALL ZGEMM( 'No transpose', 'Conjugate transpose', K, M, N, DCMPLX( -ONE ), A( M-K+1, 1 ), LDA, Q, LDA, DCMPLX( ONE ), R( M-K+1, N-M+1 ), LDA )
*
*     Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .
*
      ANORM = ZLANGE( '1', K, N, A( M-K+1, 1 ), LDA, RWORK )
      RESID = ZLANGE( '1', K, M, R( M-K+1, N-M+1 ), LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA )
      CALL ZHERK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = ZLANSY( '1', 'Upper', M, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
*     End of ZRQT02
*
      END
