      SUBROUTINE ZQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         ROGUE
      const              ROGUE = ( -1.0D+10, -1.0D+10 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )

      // Copy the first k columns of the factorization to the array Q

      CALL ZLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL ZLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )

      // Generate the first n columns of the matrix Q

      SRNAMT = 'ZUNGQR'
      CALL ZUNGQR( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy R(1:n,1:k)

      CALL ZLASET( 'Full', N, K, DCMPLX( ZERO ), DCMPLX( ZERO ), R, LDA )
      CALL ZLACPY( 'Upper', N, K, AF, LDA, R, LDA )

      // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

      CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, K, M, DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), R, LDA )

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, K, A, LDA, RWORK )
      RESID = ZLANGE( '1', N, K, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q'*Q

      CALL ZLASET( 'Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA )
      CALL ZHERK( 'Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, R, LDA )

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = ZLANSY( '1', 'Upper', N, R, LDA, RWORK )

      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS

      RETURN

      // End of ZQRT02

      }
