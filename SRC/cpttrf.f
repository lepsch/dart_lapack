      SUBROUTINE CPTTRF( N, D, E, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      COMPLEX            E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, I4;
      REAL               EII, EIR, F, G
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, MOD, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'CPTTRF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Compute the L*D*L**H (or U**H *D*U) factorization of A.

      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 20
         END IF
         EIR = REAL( E( I ) )
         EII = AIMAG( E( I ) )
         F = EIR / D( I )
         G = EII / D( I )
         E( I ) = CMPLX( F, G )
         D( I+1 ) = D( I+1 ) - F*EIR - G*EII
   10 CONTINUE

      DO 110 I = I4+1, N - 4, 4

         // Drop out of the loop if d(i) <= 0: the matrix is not positive
         // definite.

         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 20
         END IF

         // Solve for e(i) and d(i+1).

         EIR = REAL( E( I ) )
         EII = AIMAG( E( I ) )
         F = EIR / D( I )
         G = EII / D( I )
         E( I ) = CMPLX( F, G )
         D( I+1 ) = D( I+1 ) - F*EIR - G*EII

         IF( D( I+1 ).LE.ZERO ) THEN
            INFO = I+1
            GO TO 20
         END IF

         // Solve for e(i+1) and d(i+2).

         EIR = REAL( E( I+1 ) )
         EII = AIMAG( E( I+1 ) )
         F = EIR / D( I+1 )
         G = EII / D( I+1 )
         E( I+1 ) = CMPLX( F, G )
         D( I+2 ) = D( I+2 ) - F*EIR - G*EII

         IF( D( I+2 ).LE.ZERO ) THEN
            INFO = I+2
            GO TO 20
         END IF

         // Solve for e(i+2) and d(i+3).

         EIR = REAL( E( I+2 ) )
         EII = AIMAG( E( I+2 ) )
         F = EIR / D( I+2 )
         G = EII / D( I+2 )
         E( I+2 ) = CMPLX( F, G )
         D( I+3 ) = D( I+3 ) - F*EIR - G*EII

         IF( D( I+3 ).LE.ZERO ) THEN
            INFO = I+3
            GO TO 20
         END IF

         // Solve for e(i+3) and d(i+4).

         EIR = REAL( E( I+3 ) )
         EII = AIMAG( E( I+3 ) )
         F = EIR / D( I+3 )
         G = EII / D( I+3 )
         E( I+3 ) = CMPLX( F, G )
         D( I+4 ) = D( I+4 ) - F*EIR - G*EII
  110 CONTINUE

      // Check d(n) for positive definiteness.

      IF( D( N ).LE.ZERO ) INFO = N

   20 CONTINUE
      RETURN

      // End of CPTTRF

      END
