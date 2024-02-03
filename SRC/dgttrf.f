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
      if ( N.LT.0 ) {
         INFO = -1
         xerbla('DGTTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize IPIV(i) = i and DU2(I) = 0

      for (I = 1; I <= N; I++) { // 10
         IPIV( I ) = I
      } // 10
      for (I = 1; I <= N - 2; I++) { // 20
         DU2( I ) = ZERO
      } // 20

      for (I = 1; I <= N - 2; I++) { // 30
         if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {

            // No row interchange required, eliminate DL(I)

            if ( D( I ).NE.ZERO ) {
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            }
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
         }
      } // 30
      if ( N.GT.1 ) {
         I = N - 1
         if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {
            if ( D( I ).NE.ZERO ) {
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            }
         } else {
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            IPIV( I ) = I + 1
         }
      }

      // Check for a zero on the diagonal of U.

      for (I = 1; I <= N; I++) { // 40
         if ( D( I ).EQ.ZERO ) {
            INFO = I
            GO TO 50
         }
      } // 40
      } // 50

      RETURN

      // End of DGTTRF

      }
