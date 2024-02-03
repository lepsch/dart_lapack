      SUBROUTINE SGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      REAL               RESULT
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0, ONE = 1.0E0, TEN = 10.0E0 ;
      // ..
      // .. Local Scalars ..
      int                JCOL, JDIAG, JROW;
      REAL               ANORM, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY
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

         ANORM = MAX( SLANGE( '1', N, N, A, LDA, WORK ), UNFL )

         if ( ITYPE.EQ.1 ) {

            // ITYPE=1: Compute W = A - UBV'

            CALL SLACPY( ' ', N, N, A, LDA, WORK, N )
            CALL SGEMM( 'N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO, WORK( N**2+1 ), N )

            CALL SGEMM( 'N', 'C', N, N, N, -ONE, WORK( N**2+1 ), N, V, LDV, ONE, WORK, N )

         } else {

            // ITYPE=2: Compute W = A - B

            CALL SLACPY( ' ', N, N, B, LDB, WORK, N )

            DO 20 JCOL = 1, N
               DO 10 JROW = 1, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )

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

         // ITYPE=3: Compute  UU' - I

         CALL SGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N )

         DO 30 JDIAG = 1, N
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - ONE
   30    CONTINUE

         RESULT = MIN( SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), REAL( N ) ) / ( N*ULP )
      }

      RETURN

      // End of SGET51

      }
