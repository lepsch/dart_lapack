      SUBROUTINE DBDT05( M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

* ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL LSAME, IDAMAX, DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO
      IF( MIN( M, N ).LE.0 .OR. NS.LE.0 ) RETURN

      EPS = DLAMCH( 'Precision' )
      ANORM = DLANGE( 'M', M, N, A, LDA, WORK )

      // Compute U' * A * V.

      dgemm('N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO, WORK( 1+NS*NS ), M )       CALL DGEMM( 'T', 'N', NS, NS, M, -ONE, U, LDU, WORK( 1+NS*NS ), M, ZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0
      for (I = 1; I <= NS; I++) { // 10
         WORK( J+I ) =  WORK( J+I ) + S( I )
         RESID = MAX( RESID, DASUM( NS, WORK( J+1 ), 1 ) )
         J = J + NS
      } // 10

      if ( ANORM.LE.ZERO ) {
         if (RESID.NE.ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM.GE.RESID ) {
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         } else {
            if ( ANORM.LT.ONE ) {
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS )
            } else {
               RESID = MIN( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            }
         }
      }

      RETURN

      // End of DBDT05

      }
