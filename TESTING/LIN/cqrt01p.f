      SUBROUTINE CQRT01P( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            ROGUE
      const              ROGUE = ( -1.0E+10, -1.0E+10 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO, MINMN;
      REAL               ANORM, EPS, RESID
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGEQRFP, CHERK, CLACPY, CLASET, CUNGQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )

      // Factorize the matrix A in the array AF.

      SRNAMT = 'CGEQRFP'
      CALL CGEQRFP( M, N, AF, LDA, TAU, WORK, LWORK, INFO )

      // Copy details of Q

      CALL CLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      CALL CLACPY( 'Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )

      // Generate the m-by-m matrix Q

      SRNAMT = 'CUNGQR'
      CALL CUNGQR( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy R

      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA )
      CALL CLACPY( 'Upper', M, N, AF, LDA, R, LDA )

      // Compute R - Q'*A

      CALL CGEMM( 'Conjugate transpose', 'No transpose', M, N, M, CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA )

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q'*Q

      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), R, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, R, LDA )

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = CLANSY( '1', 'Upper', M, R, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS

      RETURN

      // End of CQRT01P

      }
