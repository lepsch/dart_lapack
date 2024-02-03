      SUBROUTINE CHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, LDA, LDH, LDQ, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), H( LDH, * ), Q( LDQ, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                LDWORK;
      REAL               ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CUNT01
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      END IF

      UNFL = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      OVFL = ONE / UNFL
      SMLNUM = UNFL*N / EPS

      // Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )

      // Copy A to WORK

      LDWORK = MAX( 1, N )
      CALL CLACPY( ' ', N, N, A, LDA, WORK, LDWORK )

      // Compute Q*H

      CALL CGEMM( 'No transpose', 'No transpose', N, N, N, CMPLX( ONE ), Q, LDQ, H, LDH, CMPLX( ZERO ), WORK( LDWORK*N+1 ), LDWORK )

      // Compute A - Q*H*Q'

      CALL CGEMM( 'No transpose', 'Conjugate transpose', N, N, N, CMPLX( -ONE ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, CMPLX( ONE ), WORK, LDWORK )

      ANORM = MAX( CLANGE( '1', N, N, A, LDA, RWORK ), UNFL )
      WNORM = CLANGE( '1', N, N, WORK, LDWORK, RWORK )

      // Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)

      RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N

      // Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )

      CALL CUNT01( 'Columns', N, N, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) )

      RETURN

      // End of CHST01

      }
