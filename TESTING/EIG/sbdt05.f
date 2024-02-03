      SUBROUTINE SBDT05( M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SLAMCH, SLANGE;
      // EXTERNAL LSAME, ISAMAX, SASUM, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = SLAMCH( 'Precision' );
      ANORM = SLANGE( 'M', M, N, A, LDA, WORK );

      // Compute U' * A * V.

      sgemm('N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO, WORK( 1+NS*NS ), M )       CALL SGEMM( 'T', 'N', NS, NS, M, -ONE, U, LDU, WORK( 1+NS*NS ), M, ZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK( J+I ) =  WORK( J+I ) + S( I );
         RESID = max( RESID, SASUM( NS, WORK( J+1 ), 1 ) );
         J = J + NS;
      } // 10

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS );
         } else {
            if ( ANORM < ONE ) {
               RESID = ( min( RESID, REAL( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS );
            } else {
               RESID = min( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      return;

      // End of SBDT05

      }
