      SUBROUTINE ZQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
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
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' );

      // Copy the first k columns of the factorization to the array Q

      zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      zlacpy('Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the first n columns of the matrix Q

      SRNAMT = 'ZUNGQR';
      zungqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R(1:n,1:k)

      zlaset('Full', N, K, DCMPLX( ZERO ), DCMPLX( ZERO ), R, LDA );
      zlacpy('Upper', N, K, AF, LDA, R, LDA );

      // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

      zgemm('Conjugate transpose', 'No transpose', N, K, M, DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, K, A, LDA, RWORK );
      RESID = ZLANGE( '1', N, K, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute I - Q'*Q

      zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA );
      zherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = ZLANSY( '1', 'Upper', N, R, LDA, RWORK );

      RESULT( 2 ) = ( RESID / DBLE( max( 1, M ) ) ) / EPS;

      return;

      // End of ZQRT02

      }
