      SUBROUTINE SPTTRF( N, D, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, I4;
      REAL               EI
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'SPTTRF', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Compute the L*D*L**T (or U**T*D*U) factorization of A.
*
      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
*
      DO 20 I = I4 + 1, N - 4, 4
*
         // Drop out of the loop if d(i) <= 0: the matrix is not positive
         // definite.
*
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
*
         // Solve for e(i) and d(i+1).
*
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
*
         IF( D( I+1 ).LE.ZERO ) THEN
            INFO = I + 1
            GO TO 30
         END IF
*
         // Solve for e(i+1) and d(i+2).
*
         EI = E( I+1 )
         E( I+1 ) = EI / D( I+1 )
         D( I+2 ) = D( I+2 ) - E( I+1 )*EI
*
         IF( D( I+2 ).LE.ZERO ) THEN
            INFO = I + 2
            GO TO 30
         END IF
*
         // Solve for e(i+2) and d(i+3).
*
         EI = E( I+2 )
         E( I+2 ) = EI / D( I+2 )
         D( I+3 ) = D( I+3 ) - E( I+2 )*EI
*
         IF( D( I+3 ).LE.ZERO ) THEN
            INFO = I + 3
            GO TO 30
         END IF
*
         // Solve for e(i+3) and d(i+4).
*
         EI = E( I+3 )
         E( I+3 ) = EI / D( I+3 )
         D( I+4 ) = D( I+4 ) - E( I+3 )*EI
   20 CONTINUE
*
      // Check d(n) for positive definiteness.
*
      IF( D( N ).LE.ZERO ) INFO = N
*
   30 CONTINUE
      RETURN
*
      // End of SPTTRF
*
      END
