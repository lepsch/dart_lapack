      SUBROUTINE CQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
*     ..
*     .. Local Scalars ..
      int                INFO
      REAL               ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      EXTERNAL           CLANGE, CLANSY, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CHERK, CLACPY, CLASET, CUNGQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, REAL
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the first k columns of the factorization to the array Q
*
      CALL CLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL CLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the first n columns of the matrix Q
*
      SRNAMT = 'CUNGQR'
      CALL CUNGQR( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R(1:n,1:k)
*
      CALL CLASET( 'Full', N, K, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA )
      CALL CLACPY( 'Upper', N, K, AF, LDA, R, LDA )
*
*     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', N, K, M, CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', M, K, A, LDA, RWORK )
      RESID = CLANGE( '1', N, K, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of CQRT02
*
      END
