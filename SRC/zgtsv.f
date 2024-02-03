      SUBROUTINE ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
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
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGTSV ', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 ) RETURN
*
      DO 30 K = 1, N - 1
         IF( DL( K ).EQ.ZERO ) THEN
*
            // Subdiagonal is zero, no elimination is required.
*
            IF( D( K ).EQ.ZERO ) THEN
*
               // Diagonal is zero: set INFO = K and return; a unique
               // solution can not be found.
*
               INFO = K
               RETURN
            END IF
         ELSE IF( CABS1( D( K ) ).GE.CABS1( DL( K ) ) ) THEN
*
            // No row interchange required
*
            MULT = DL( K ) / D( K )
            D( K+1 ) = D( K+1 ) - MULT*DU( K )
            DO 10 J = 1, NRHS
               B( K+1, J ) = B( K+1, J ) - MULT*B( K, J )
   10       CONTINUE
            IF( K.LT.( N-1 ) ) DL( K ) = ZERO
         ELSE
*
            // Interchange rows K and K+1
*
            MULT = D( K ) / DL( K )
            D( K ) = DL( K )
            TEMP = D( K+1 )
            D( K+1 ) = DU( K ) - MULT*TEMP
            IF( K.LT.( N-1 ) ) THEN
               DL( K ) = DU( K+1 )
               DU( K+1 ) = -MULT*DL( K )
            END IF
            DU( K ) = TEMP
            DO 20 J = 1, NRHS
               TEMP = B( K, J )
               B( K, J ) = B( K+1, J )
               B( K+1, J ) = TEMP - MULT*B( K+1, J )
   20       CONTINUE
         END IF
   30 CONTINUE
      IF( D( N ).EQ.ZERO ) THEN
         INFO = N
         RETURN
      END IF
*
      // Back solve with the matrix U from the factorization.
*
      DO 50 J = 1, NRHS
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 40 K = N - 2, 1, -1
            B( K, J ) = ( B( K, J )-DU( K )*B( K+1, J )-DL( K )* B( K+2, J ) ) / D( K )
   40    CONTINUE
   50 CONTINUE
*
      RETURN
*
      // End of ZGTSV
*
      END
