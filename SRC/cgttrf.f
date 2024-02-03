      SUBROUTINE CGTTRF( N, DL, D, DU, DU2, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            D( * ), DL( * ), DU( * ), DU2( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      COMPLEX            FACT, TEMP, ZDUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
         xerbla('CGTTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize IPIV(i) = i and DU2(i) = 0

      DO 10 I = 1, N
         IPIV( I ) = I
   10 CONTINUE
      DO 20 I = 1, N - 2
         DU2( I ) = ZERO
   20 CONTINUE

      DO 30 I = 1, N - 2
         if ( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) {

            // No row interchange required, eliminate DL(I)

            if ( CABS1( D( I ) ).NE.ZERO ) {
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
   30 CONTINUE
      if ( N.GT.1 ) {
         I = N - 1
         if ( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) {
            if ( CABS1( D( I ) ).NE.ZERO ) {
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

      DO 40 I = 1, N
         if ( CABS1( D( I ) ).EQ.ZERO ) {
            INFO = I
            GO TO 50
         }
   40 CONTINUE
   50 CONTINUE

      RETURN

      // End of CGTTRF

      }
