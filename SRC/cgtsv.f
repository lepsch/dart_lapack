      SUBROUTINE CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, K;
      COMPLEX            MULT, TEMP, ZDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
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
      } else if ( NRHS.LT.0 ) {
         INFO = -2
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         xerbla('CGTSV ', -INFO );
         RETURN
      }

      if (N == 0) RETURN;

      for (K = 1; K <= N - 1; K++) { // 30
         if ( DL( K ) == ZERO ) {

            // Subdiagonal is zero, no elimination is required.

            if ( D( K ) == ZERO ) {

               // Diagonal is zero: set INFO = K and return; a unique
               // solution can not be found.

               INFO = K
               RETURN
            }
         } else if ( CABS1( D( K ) ).GE.CABS1( DL( K ) ) ) {

            // No row interchange required

            MULT = DL( K ) / D( K )
            D( K+1 ) = D( K+1 ) - MULT*DU( K )
            for (J = 1; J <= NRHS; J++) { // 10
               B( K+1, J ) = B( K+1, J ) - MULT*B( K, J )
            } // 10
            IF( K.LT.( N-1 ) ) DL( K ) = ZERO
         } else {

            // Interchange rows K and K+1

            MULT = D( K ) / DL( K )
            D( K ) = DL( K )
            TEMP = D( K+1 )
            D( K+1 ) = DU( K ) - MULT*TEMP
            if ( K.LT.( N-1 ) ) {
               DL( K ) = DU( K+1 )
               DU( K+1 ) = -MULT*DL( K )
            }
            DU( K ) = TEMP
            for (J = 1; J <= NRHS; J++) { // 20
               TEMP = B( K, J )
               B( K, J ) = B( K+1, J )
               B( K+1, J ) = TEMP - MULT*B( K+1, J )
            } // 20
         }
      } // 30
      if ( D( N ) == ZERO ) {
         INFO = N
         RETURN
      }

      // Back solve with the matrix U from the factorization.

      for (J = 1; J <= NRHS; J++) { // 50
         B( N, J ) = B( N, J ) / D( N )
         if (N.GT.1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
         DO 40 K = N - 2, 1, -1
            B( K, J ) = ( B( K, J )-DU( K )*B( K+1, J )-DL( K )* B( K+2, J ) ) / D( K )
         } // 40
      } // 50

      RETURN

      // End of CGTSV

      }
