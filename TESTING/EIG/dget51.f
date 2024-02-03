      SUBROUTINE DGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 ;
      // ..
      // .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      double             ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT = ZERO
      if (N.LE.0) RETURN;

      // Constants

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // Some Error Checks

      if ( ITYPE < 1 || ITYPE.GT.3 ) {
         RESULT = TEN / ULP
         RETURN
      }

      if ( ITYPE.LE.2 ) {

         // Tests scaled by the norm(A)

         ANORM = MAX( DLANGE( '1', N, N, A, LDA, WORK ), UNFL )

         if ( ITYPE == 1 ) {

            // ITYPE=1: Compute W = A - UBV'

            dlacpy(' ', N, N, A, LDA, WORK, N );
            dgemm('N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO, WORK( N**2+1 ), N );

            dgemm('N', 'C', N, N, N, -ONE, WORK( N**2+1 ), N, V, LDV, ONE, WORK, N );

         } else {

            // ITYPE=2: Compute W = A - B

            dlacpy(' ', N, N, B, LDB, WORK, N );

            for (JCOL = 1; JCOL <= N; JCOL++) { // 20
               for (JROW = 1; JROW <= N; JROW++) { // 10
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
               } // 10
            } // 20
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )

         if ( ANORM.GT.WNORM ) {
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         } else {
            if ( ANORM < ONE ) {
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            } else {
               RESULT = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            }
         }

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  UU' - I

         dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

         for (JDIAG = 1; JDIAG <= N; JDIAG++) { // 30
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - ONE
         } // 30

         RESULT = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP )
      }

      RETURN

      // End of DGET51

      }
