      SUBROUTINE SGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      REAL               RESULT;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      // ..
      // .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      REAL               ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT = ZERO;
      if (N <= 0) RETURN;

      // Constants

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      // Some Error Checks

      if ( ITYPE < 1 || ITYPE > 3 ) {
         RESULT = TEN / ULP;
         return;
      }

      if ( ITYPE <= 2 ) {

         // Tests scaled by the norm(A)

         ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), UNFL );

         if ( ITYPE == 1 ) {

            // ITYPE=1: Compute W = A - UBV'

            slacpy(' ', N, N, A, LDA, WORK, N );
            sgemm('N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO, WORK( N**2+1 ), N );

            sgemm('N', 'C', N, N, N, -ONE, WORK( N**2+1 ), N, V, LDV, ONE, WORK, N );

         } else {

            // ITYPE=2: Compute W = A - B

            slacpy(' ', N, N, B, LDB, WORK, N );

            for (JCOL = 1; JCOL <= N; JCOL++) { // 20
               for (JROW = 1; JROW <= N; JROW++) { // 10
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL );
               } // 10
            } // 20
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) );

         if ( ANORM > WNORM ) {
            RESULT = ( WNORM / ANORM ) / ( N*ULP );
         } else {
            if ( ANORM < ONE ) {
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
            } else {
               RESULT = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP );
            }
         }

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  UU' - I

         sgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

         for (JDIAG = 1; JDIAG <= N; JDIAG++) { // 30
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - ONE;
         } // 30

         RESULT = MIN( SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP );
      }

      return;

      // End of SGET51

      }
