      void sgttrf(N, DL, D, DU, DU2, IPIV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double               D( * ), DL( * ), DU( * ), DU2( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double               FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('SGTTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize IPIV(i) = i and DU2(I) = 0

      for (I = 1; I <= N; I++) { // 10
         IPIV[I] = I;
      } // 10
      for (I = 1; I <= N - 2; I++) { // 20
         DU2[I] = ZERO;
      } // 20

      for (I = 1; I <= N - 2; I++) { // 30
         if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {

            // No row interchange required, eliminate DL(I)

            if ( D( I ) != ZERO ) {
               FACT = DL( I ) / D( I );
               DL[I] = FACT;
               D[I+1] = D( I+1 ) - FACT*DU( I );
            }
         } else {

            // Interchange rows I and I+1, eliminate DL(I)

            FACT = D( I ) / DL( I );
            D[I] = DL( I );
            DL[I] = FACT;
            TEMP = DU( I );
            DU[I] = D( I+1 );
            D[I+1] = TEMP - FACT*D( I+1 );
            DU2[I] = DU( I+1 );
            DU[I+1] = -FACT*DU( I+1 );
            IPIV[I] = I + 1;
         }
      } // 30
      if ( N > 1 ) {
         I = N - 1;
         if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {
            if ( D( I ) != ZERO ) {
               FACT = DL( I ) / D( I );
               DL[I] = FACT;
               D[I+1] = D( I+1 ) - FACT*DU( I );
            }
         } else {
            FACT = D( I ) / DL( I );
            D[I] = DL( I );
            DL[I] = FACT;
            TEMP = DU( I );
            DU[I] = D( I+1 );
            D[I+1] = TEMP - FACT*D( I+1 );
            IPIV[I] = I + 1;
         }
      }

      // Check for a zero on the diagonal of U.

      for (I = 1; I <= N; I++) { // 40
         if ( D( I ) == ZERO ) {
            INFO = I;
            GO TO 50;
         }
      } // 40
      } // 50

      return;
      }
