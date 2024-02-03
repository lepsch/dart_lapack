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
      IF( N.LE.0 ) RETURN

      // Constants

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // Some Error Checks

      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         RESULT = TEN / ULP
         RETURN
      END IF

      IF( ITYPE.LE.2 ) THEN

         // Tests scaled by the norm(A)

         ANORM = MAX( DLANGE( '1', N, N, A, LDA, WORK ), UNFL )

         IF( ITYPE.EQ.1 ) THEN

            // ITYPE=1: Compute W = A - UBV'

            CALL DLACPY( ' ', N, N, A, LDA, WORK, N )
            CALL DGEMM( 'N', 'N', N, N, N, ONE, U, LDU, B, LDB, ZERO, WORK( N**2+1 ), N )

            CALL DGEMM( 'N', 'C', N, N, N, -ONE, WORK( N**2+1 ), N, V, LDV, ONE, WORK, N )

         } else {

            // ITYPE=2: Compute W = A - B

            CALL DLACPY( ' ', N, N, B, LDB, WORK, N )

            DO 20 JCOL = 1, N
               DO 10 JROW = 1, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         END IF

         // Compute norm(W)/ ( ulp*norm(A) )

         WNORM = DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )

         IF( ANORM.GT.WNORM ) THEN
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         } else {
            IF( ANORM.LT.ONE ) THEN
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            } else {
               RESULT = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            END IF
         END IF

      } else {

         // Tests not scaled by norm(A)

         // ITYPE=3: Compute  UU' - I

         CALL DGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N )

         DO 30 JDIAG = 1, N
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+ 1 ) - ONE
   30    CONTINUE

         RESULT = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP )
      END IF

      RETURN

      // End of DGET51

      }
