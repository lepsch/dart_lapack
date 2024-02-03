      SUBROUTINE CBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDC, LDU, M, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * )
      // ..

* ======================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               BNORM, EPS, REALMN
      // ..
      // .. External Functions ..
      REAL               CLANGE, SCASUM, SLAMCH
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
      REALMN = REAL( MAX( M, N ) )
      EPS = SLAMCH( 'Precision' )

      // Compute norm(B - U * C)

      DO 10 J = 1, N
         CALL CCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL CGEMV( 'No transpose', M, M, -CMPLX( ONE ), U, LDU, C( 1, J ), 1, CMPLX( ONE ), WORK, 1 )
         RESID = MAX( RESID, SCASUM( M, WORK, 1 ) )
   10 CONTINUE

      // Compute norm of B.

      BNORM = CLANGE( '1', M, N, B, LDB, RWORK )

      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REALMN*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) / ( REALMN*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
            END IF
         END IF
      END IF
      RETURN

      // End of CBDT02

      END
