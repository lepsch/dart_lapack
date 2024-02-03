      SUBROUTINE DBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDB, LDC, LDU, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
      // ..
*
* ======================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                J;
      double             BNORM, EPS, REALMN;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Quick return if possible
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
      REALMN = DBLE( MAX( M, N ) )
      EPS = DLAMCH( 'Precision' )
*
      // Compute norm(B - U * C)
*
      DO 10 J = 1, N
         CALL DCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL DGEMV( 'No transpose', M, M, -ONE, U, LDU, C( 1, J ), 1, ONE, WORK, 1 )
         RESID = MAX( RESID, DASUM( M, WORK, 1 ) )
   10 CONTINUE
*
      // Compute norm of B.
*
      BNORM = DLANGE( '1', M, N, B, LDB, WORK )
*
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
*
      // End of DBDT02
*
      END
