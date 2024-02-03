      SUBROUTINE ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         D( * ), DL( * ), DU( * ), DU2( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX*16         FACT, TEMP, ZDUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N < 0 ) {
         INFO = -1
         xerbla('ZGTTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Initialize IPIV(i) = i and DU2(i) = 0

      for (I = 1; I <= N; I++) { // 10
         IPIV( I ) = I
      } // 10
      for (I = 1; I <= N - 2; I++) { // 20
         DU2( I ) = ZERO
      } // 20

      for (I = 1; I <= N - 2; I++) { // 30
         if ( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) {

            // No row interchange required, eliminate DL(I)

            if ( CABS1( D( I ) ) != ZERO ) {
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
      if ( N > 1 ) {
         I = N - 1
         if ( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) {
            if ( CABS1( D( I ) ) != ZERO ) {
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
         if ( CABS1( D( I ) ) == ZERO ) {
            INFO = I
            GO TO 50
         }
      } // 40
      } // 50

      RETURN

      // End of ZGTTRF

      }
