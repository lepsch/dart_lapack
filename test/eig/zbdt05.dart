      void zbdt05(final int M, final int N, final Matrix<double> A, final int LDA, final int S, final int NS, final Matrix<double> U, final int LDU, final Matrix<double> VT, final int LDVT, final Array<double> _WORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDU, LDVT, M, N, NS;
      double             RESID;
      double             S( * );
      Complex         A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, J;
      double             ANORM, EPS;
      double             DUM( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DASUM, DZASUM, DLAMCH, ZLANGE;
      // EXTERNAL lsame, idamax, DASUM, DZASUM, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = dlamch( 'Precision' );
      ANORM = ZLANGE( 'M', M, N, A, LDA, DUM );

      // Compute U' * A * V.

      zgemm('N', 'C', M, NS, N, CONE, A, LDA, VT, LDVT, CZERO, WORK( 1+NS*NS ), M )       CALL ZGEMM( 'C', 'N', NS, NS, M, -CONE, U, LDU, WORK( 1+NS*NS ), M, CZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK[J+I] = WORK( J+I ) + DCMPLX( S( I ), ZERO );
         RESID = max( RESID, DZASUM( NS, WORK( J+1 ), 1 ) );
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

      }
