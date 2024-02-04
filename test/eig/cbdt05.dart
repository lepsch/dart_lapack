      void cbdt05(M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               S( * );
      COMPLEX            A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS;
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SASUM, SCASUM, SLAMCH, CLANGE;
      // EXTERNAL lsame, ISAMAX, SASUM, SCASUM, SLAMCH, CLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if( min( M, N ) <= 0 || NS <= 0 ) return;

      EPS = SLAMCH( 'Precision' );
      ANORM = CLANGE( 'M', M, N, A, LDA, DUM );

      // Compute U' * A * V.

      cgemm('N', 'C', M, NS, N, CONE, A, LDA, VT, LDVT, CZERO, WORK( 1+NS*NS ), M )       CALL CGEMM( 'C', 'N', NS, NS, M, -CONE, U, LDU, WORK( 1+NS*NS ), M, CZERO, WORK, NS );

      // norm(S - U' * B * V)

      J = 0;
      for (I = 1; I <= NS; I++) { // 10
         WORK[J+I] = WORK( J+I ) + CMPLX( S( I ), ZERO );
         RESID = max( RESID, SCASUM( NS, WORK( J+1 ), 1 ) );
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
      }
