      SUBROUTINE ZGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITYPE, LDA, LDB, LDU, LDV, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 10.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
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

      RESULT = ZERO
      IF( N.LE.0 ) RETURN

      // Constants

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // Some Error Checks

      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         RESULT = TEN / ULP
         RETURN
      }

      if ( ITYPE.LE.2 ) {

         // Tests scaled by the norm(A)

         ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL )

         if ( ITYPE.EQ.1 ) {

            // ITYPE=1: Compute W = A - U B V**H

            zlacpy(' ', N, N, A, LDA, WORK, N );
            zgemm('N', 'N', N, N, N, CONE, U, LDU, B, LDB, CZERO, WORK( N**2+1 ), N );

            zgemm('N', 'C', N, N, N, -CONE, WORK( N**2+1 ), N, V, LDV, CONE, WORK, N );

         } else {

            // ITYPE=2: Compute W = A - B

            zlacpy(' ', N, N, B, LDB, WORK, N );

            DO 20 JCOL = 1, N
               DO 10 JROW = 1, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         }

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = ZLANGE( '1', N, N, WORK, N, RWORK )

         if ( ANORM.GT.WNORM ) {
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         } else {
            if ( ANORM.LT.ONE ) {
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            } else {
               RESULT = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            }
         }

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  U U**H - I

         zgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

         DO 30 JDIAG = 1, N
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - CONE
   30    CONTINUE

         RESULT = MIN( ZLANGE( '1', N, N, WORK, N, RWORK ), DBLE( N ) ) / ( N*ULP )
      }

      RETURN

      // End of ZGET51

      }
