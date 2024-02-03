      void zbdt05(M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             S( * );
      Complex         A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ANORM, EPS;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DASUM, DZASUM, DLAMCH, ZLANGE;
      // EXTERNAL LSAME, IDAMAX, DASUM, DZASUM, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = DLAMCH( 'Precision' );
      ANORM = ZLANGE( 'M', M, N, A, LDA, DUM );

      // Compute U' * A * V.

      zgemm('N', 'C', M, NS, N, CONE, A, LDA, VT, LDVT, CZERO, WORK( 1+NS*NS ), M )       CALL ZGEMM( 'C', 'N', NS, NS, M, -CONE, U, LDU, WORK( 1+NS*NS ), M, CZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK( J+I ) =  WORK( J+I ) + DCMPLX( S( I ), ZERO );
         RESID = max( RESID, DZASUM( NS, WORK( J+1 ), 1 ) );
         J = J + NS;
      } // 10

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( ANORM >= RESID ) {
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS );
         } else {
            if ( ANORM < ONE ) {
               RESID = ( min( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS );
            } else {
               RESID = min( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS );
            }
         }
      }

      return;

      // End of ZBDT05

      }
