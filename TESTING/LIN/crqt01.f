      SUBROUTINE CRQT01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )

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
      // EXTERNAL CGEMM, CGERQF, CHERK, CLACPY, CLASET, CUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

      SRNAMT = 'CGERQF'
      cgerqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      claset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M.LE.N ) {
         if (M.GT.0 .AND. M.LT.N) CALL CLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL CLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N.GT.1) CALL CLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }

      // Generate the n-by-n matrix Q

      SRNAMT = 'CUNGRQ'
      cungrq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA );
      if ( M.LE.N ) {
         if (M.GT.0) CALL CLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         if (M.GT.N .AND. N.GT.0) CALL CLACPY( 'Full', M-N, N, AF, LDA, R, LDA )          IF( N.GT.0 ) CALL CLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Compute R - A*Q'

      cgemm('No transpose', 'Conjugate transpose', M, N, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), R, LDA );

      // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, R, LDA, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute I - Q*Q'

      claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA );
      cherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of CRQT01

      }
