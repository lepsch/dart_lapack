      void zget51(ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      double             ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT = ZERO;
      if (N <= 0) return;

      // Constants

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );

      // Some Error Checks

      if ( ITYPE < 1 || ITYPE > 3 ) {
         RESULT = TEN / ULP;
         return;
      }

      if ( ITYPE <= 2 ) {

         // Tests scaled by the norm(A)

         ANORM = max( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL );

         if ( ITYPE == 1 ) {

            // ITYPE=1: Compute W = A - U B V**H

            zlacpy(' ', N, N, A, LDA, WORK, N );
            zgemm('N', 'N', N, N, N, CONE, U, LDU, B, LDB, CZERO, WORK( N**2+1 ), N );

            zgemm('N', 'C', N, N, N, -CONE, WORK( N**2+1 ), N, V, LDV, CONE, WORK, N );

         } else {

            // ITYPE=2: Compute W = A - B

            zlacpy(' ', N, N, B, LDB, WORK, N );

            for (JCOL = 1; JCOL <= N; JCOL++) { // 20
               for (JROW = 1; JROW <= N; JROW++) { // 10
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL );
               } // 10
            } // 20
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = ZLANGE( '1', N, N, WORK, N, RWORK );

         if ( ANORM > WNORM ) {
            RESULT = ( WNORM / ANORM ) / ( N*ULP );
         } else {
            if ( ANORM < ONE ) {
               RESULT = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
            } else {
               RESULT = min( WNORM / ANORM, DBLE( N ) ) / ( N*ULP );
            }
         }

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  U U**H - I

         zgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

         for (JDIAG = 1; JDIAG <= N; JDIAG++) { // 30
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - CONE;
         } // 30

         RESULT = min( ZLANGE( '1', N, N, WORK, N, RWORK ), DBLE( N ) ) / ( N*ULP );
      }

      return;
      }
