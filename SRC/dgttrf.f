      SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             D( * ), DL( * ), DU( * ), DU2( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DGTTRF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize IPIV(i) = i and DU2(I) = 0

      DO 10 I = 1, N
         IPIV( I ) = I
   10 CONTINUE
      DO 20 I = 1, N - 2
         DU2( I ) = ZERO
   20 CONTINUE

      DO 30 I = 1, N - 2
         IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN

            // No row interchange required, eliminate DL(I)

            IF( D( I ).NE.ZERO ) THEN
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            END IF
         } else {

            // Interchange rows I and I+1, eliminate DL(I)

            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            DU2( I ) = DU( I+1 )
            DU( I+1 ) = -FACT*DU( I+1 )
            IPIV( I ) = I + 1
         END IF
   30 CONTINUE
      IF( N.GT.1 ) THEN
         I = N - 1
         IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
            IF( D( I ).NE.ZERO ) THEN
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            END IF
         } else {
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            IPIV( I ) = I + 1
         END IF
      END IF

      // Check for a zero on the diagonal of U.

      DO 40 I = 1, N
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I
            GO TO 50
         END IF
   40 CONTINUE
   50 CONTINUE

      RETURN

      // End of DGTTRF

      }
