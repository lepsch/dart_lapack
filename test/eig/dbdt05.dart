      void dbdt05(M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL lsame, idamax, DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = dlamch( 'Precision' );
      ANORM = dlange( 'M', M, N, A, LDA, WORK );

      // Compute U' * A * V.

      dgemm('N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO, WORK( 1+NS*NS ), M )       CALL DGEMM( 'T', 'N', NS, NS, M, -ONE, U, LDU, WORK( 1+NS*NS ), M, ZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK[J+I] = WORK( J+I ) + S( I );
         RESID = max( RESID, dasum( NS, WORK( J+1 ), 1 ) );
         J = J + NS;
      } // 10

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( N.toDouble()*EPS );
         } else {
            if ( ANORM < ONE ) {
               RESID = ( min( RESID, (N).toDouble()*ANORM ) / ANORM ) / ( N.toDouble()*EPS );
            } else {
               RESID = min( RESID / ANORM, (N).toDouble() ) / ( N.toDouble()*EPS );
            }
         }
      }

      return;
      }
