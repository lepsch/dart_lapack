      SUBROUTINE CGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      REAL               RESULT
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      REAL               ANORM, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT = ZERO
      IF( N.LE.0 ) RETURN

      // Constants

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      // Some Error Checks

      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         RESULT = TEN / ULP
         RETURN
      }

      if ( ITYPE.LE.2 ) {

         // Tests scaled by the norm(A)

         ANORM = MAX( CLANGE( '1', N, N, A, LDA, RWORK ), UNFL )

         if ( ITYPE.EQ.1 ) {

            // ITYPE=1: Compute W = A - U B V**H

            clacpy(' ', N, N, A, LDA, WORK, N );
            cgemm('N', 'N', N, N, N, CONE, U, LDU, B, LDB, CZERO, WORK( N**2+1 ), N );

            cgemm('N', 'C', N, N, N, -CONE, WORK( N**2+1 ), N, V, LDV, CONE, WORK, N );

         } else {

            // ITYPE=2: Compute W = A - B

            clacpy(' ', N, N, B, LDB, WORK, N );

            for (JCOL = 1; JCOL <= N; JCOL++) { // 20
               for (JROW = 1; JROW <= N; JROW++) { // 10
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = CLANGE( '1', N, N, WORK, N, RWORK )

         if ( ANORM.GT.WNORM ) {
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         } else {
            if ( ANORM.LT.ONE ) {
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            } else {
               RESULT = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
            }
         }

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  U U**H - I

         cgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

         for (JDIAG = 1; JDIAG <= N; JDIAG++) { // 30
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - CONE
   30    CONTINUE

         RESULT = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP )
      }

      RETURN

      // End of CGET51

      }
