      SUBROUTINE ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, K;
      COMPLEX*16         MULT, TEMP, ZDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
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
         xerbla('ZGTSV ', -INFO );
         RETURN
      }

      IF( N.EQ.0 ) RETURN

      DO 30 K = 1, N - 1
         if ( DL( K ).EQ.ZERO ) {

            // Subdiagonal is zero, no elimination is required.

            if ( D( K ).EQ.ZERO ) {

               // Diagonal is zero: set INFO = K and return; a unique
               // solution can not be found.

               INFO = K
               RETURN
            }
         } else if ( CABS1( D( K ) ).GE.CABS1( DL( K ) ) ) {

            // No row interchange required

            MULT = DL( K ) / D( K )
            D( K+1 ) = D( K+1 ) - MULT*DU( K )
            DO 10 J = 1, NRHS
               B( K+1, J ) = B( K+1, J ) - MULT*B( K, J )
   10       CONTINUE
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
            DO 20 J = 1, NRHS
               TEMP = B( K, J )
               B( K, J ) = B( K+1, J )
               B( K+1, J ) = TEMP - MULT*B( K+1, J )
   20       CONTINUE
         }
   30 CONTINUE
      if ( D( N ).EQ.ZERO ) {
         INFO = N
         RETURN
      }

      // Back solve with the matrix U from the factorization.

      DO 50 J = 1, NRHS
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 40 K = N - 2, 1, -1
            B( K, J ) = ( B( K, J )-DU( K )*B( K+1, J )-DL( K )* B( K+2, J ) ) / D( K )
   40    CONTINUE
   50 CONTINUE

      RETURN

      // End of ZGTSV

      }
