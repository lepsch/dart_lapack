      SUBROUTINE ZQRT01P( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         ROGUE
      PARAMETER          ( ROGUE = ( -1.0D+10, -1.0D+10 ) )
*     ..
*     .. Local Scalars ..
      int                INFO, MINMN;
      double             ANORM, EPS, RESID;
*     ..
*     .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANSY;
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZGEQRFP, ZHERK, ZLACPY, ZLASET, ZUNGQR
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Scalars in Common ..
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'ZGEQRFP'
      CALL ZGEQRFP( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL ZLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      CALL ZLACPY( 'Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'ZUNGQR'
      CALL ZUNGQR( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), R, LDA )
      CALL ZLACPY( 'Upper', M, N, AF, LDA, R, LDA )
*
*     Compute R - Q'*A
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M, DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), R, LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      RESID = ZLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA )
      CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = ZLANSY( '1', 'Upper', M, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of ZQRT01P
*
      END
