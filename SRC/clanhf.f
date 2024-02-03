      REAL FUNCTION CLANHF( NORM, TRANSR, UPLO, N, A, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, TRANSR, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      REAL               WORK( 0: * )
      COMPLEX            A( 0: * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, IFM, ILU, NOE, N1, K, L, LDA;
      REAL               SCALE, S, VALUE, AA, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         CLANHF = ZERO
         RETURN
      } else if ( N.EQ.1 ) {
         CLANHF = ABS(REAL(A(0)))
         RETURN
      }

      // set noe = 1 if n is odd. if n is even set noe=0

      NOE = 1
      IF( MOD( N, 2 ).EQ.0 ) NOE = 0

      // set ifm = 0 when form='C' or 'c' and 1 otherwise

      IFM = 1
      IF( LSAME( TRANSR, 'C' ) ) IFM = 0

      // set ilu = 0 when uplo='U or 'u' and 1 otherwise

      ILU = 1
      IF( LSAME( UPLO, 'U' ) ) ILU = 0

      // set lda = (n+1)/2 when ifm = 0
      // set lda = n when ifm = 1 and noe = 1
      // set lda = n+1 when ifm = 1 and noe = 0

      if ( IFM.EQ.1 ) {
         if ( NOE.EQ.1 ) {
            LDA = N
         } else {
            // noe=0
            LDA = N + 1
         }
      } else {
         // ifm=0
         LDA = ( N+1 ) / 2
      }

      if ( LSAME( NORM, 'M' ) ) {

        // Find max(abs(A(i,j))).

         K = ( N+1 ) / 2
         VALUE = ZERO
         if ( NOE.EQ.1 ) {
            // n is odd & n = k + k - 1
            if ( IFM.EQ.1 ) {
               // A is n by k
               if ( ILU.EQ.1 ) {
                  // uplo ='L'
                  J = 0
                  // -> L(0,0)
                  TEMP = ABS( REAL( A( J+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO I = 1, N - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  DO J = 1, K - 1
                     DO I = 0, J - 2
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J - 1
                     // L(k+j,k+j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J
                     // -> L(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J + 1, N - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
               } else {
                  // uplo = 'U'
                  DO J = 0, K - 2
                     DO I = 0, K + J - 2
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = K + J - 1
                     // -> U(i,i)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = I + 1
                     // =k+j; i -> U(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = K + J + 1, N - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  DO I = 0, N - 2
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     // j=k-1
                  END DO
                  // i=n-1 -> U(n-1,n-1)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
               }
            } else {
               // xpose case; A is k by n
               if ( ILU.EQ.1 ) {
                  // uplo ='L'
                  DO J = 0, K - 2
                     DO I = 0, J - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J
                     // L(i,i)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J + 1
                     // L(j+k,j+k)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J + 2, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  J = K - 1
                  DO I = 0, K - 2
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  I = K - 1
                  // -> L(i,i) is at A(i,j)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO J = K, N - 1
                     DO I = 0, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
               } else {
                  // uplo = 'U'
                  DO J = 0, K - 2
                     DO I = 0, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  J = K - 1
                  // -> U(j,j) is at A(0,j)
                  TEMP = ABS( REAL( A( 0+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO I = 1, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  DO J = K, N - 1
                     DO I = 0, J - K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J - K
                     // -> U(i,i) at A(i,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J - K + 1
                     // U(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J - K + 2, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
               }
            }
         } else {
            // n is even & k = n/2
            if ( IFM.EQ.1 ) {
               // A is n+1 by k
               if ( ILU.EQ.1 ) {
                  // uplo ='L'
                  J = 0
                  // -> L(k,k) & j=1 -> L(0,0)
                  TEMP = ABS( REAL( A( J+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  TEMP = ABS( REAL( A( J+1+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO I = 2, N
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  DO J = 1, K - 1
                     DO I = 0, J - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J
                     // L(k+j,k+j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J + 1
                     // -> L(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J + 2, N
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
               } else {
                  // uplo = 'U'
                  DO J = 0, K - 2
                     DO I = 0, K + J - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = K + J
                     // -> U(i,i)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = I + 1
                     // =k+j+1; i -> U(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = K + J + 2, N
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  DO I = 0, N - 2
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  // j=k-1
                  END DO
                  // i=n-1 -> U(n-1,n-1)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  I = N
                  // -> U(k-1,k-1)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
               }
            } else {
               // xpose case; A is k by n+1
               if ( ILU.EQ.1 ) {
                  // uplo ='L'
                  J = 0
                  // -> L(k,k) at A(0,0)
                  TEMP = ABS( REAL( A( J+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO I = 1, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  DO J = 1, K - 1
                     DO I = 0, J - 2
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J - 1
                     // L(i,i)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J
                     // L(j+k,j+k)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J + 1, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  J = K
                  DO I = 0, K - 2
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  I = K - 1
                  // -> L(i,i) is at A(i,j)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO J = K + 1, N
                     DO I = 0, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
               } else {
                  // uplo = 'U'
                  DO J = 0, K - 1
                     DO I = 0, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  J = K
                  // -> U(j,j) is at A(0,j)
                  TEMP = ABS( REAL( A( 0+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  DO I = 1, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  DO J = K + 1, N - 1
                     DO I = 0, J - K - 2
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                     I = J - K - 1
                     // -> U(i,i) at A(i,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     I = J - K
                     // U(j,j)
                     TEMP = ABS( REAL( A( I+J*LDA ) ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     DO I = J - K + 1, K - 1
                        TEMP = ABS( A( I+J*LDA ) )
                        IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                     END DO
                  END DO
                  J = N
                  DO I = 0, K - 2
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
                  I = K - 1
                  // U(k,k) at A(i,j)
                  TEMP = ABS( REAL( A( I+J*LDA ) ) )
                  IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
               }
            }
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

        // Find normI(A) ( = norm1(A), since A is Hermitian).

         if ( IFM.EQ.1 ) {
            // A is 'N'
            K = N / 2
            if ( NOE.EQ.1 ) {
               // n is odd & A is n by (n+1)/2
               if ( ILU.EQ.0 ) {
                  // uplo = 'U'
                  DO I = 0, K - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K
                     S = ZERO
                     DO I = 0, K + J - 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     END DO
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     IF( I.EQ.K+K ) GO TO 10
                     I = I + 1
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     DO L = J + 1, K - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
   10             CONTINUE
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu = 1 & uplo = 'L'
                  K = K + 1
                  // k=(n+1)/2 for n odd and ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = K - 1, 0, -1
                     S = ZERO
                     DO I = 0, J - 2
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     END DO
                     if ( J.GT.0 ) {
                        AA = ABS( REAL( A( I+J*LDA ) ) )
                        // -> A(j+k,j+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + S
                        // i=j
                        I = I + 1
                     }
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j,j)
                     WORK( J ) = AA
                     S = ZERO
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            } else {
               // n is even & A is n+1 by k = n/2
               if ( ILU.EQ.0 ) {
                  // uplo = 'U'
                  DO I = 0, K - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 1
                     S = ZERO
                     DO I = 0, K + J - 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     END DO
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     I = I + 1
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     DO L = J + 1, K - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu = 1 & uplo = 'L'
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = K - 1, 0, -1
                     S = ZERO
                     DO I = 0, J - 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     END DO
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j+k,j+k)
                     S = S + AA
                     WORK( I+K ) = WORK( I+K ) + S
                     // i=j
                     I = I + 1
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // -> A(j,j)
                     WORK( J ) = AA
                     S = ZERO
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            }
         } else {
            // ifm=0
            K = N / 2
            if ( NOE.EQ.1 ) {
               // n is odd & A is (n+1)/2 by n
               if ( ILU.EQ.0 ) {
                  // uplo = 'U'
                  N1 = K
                  // n/2
                  K = K + 1
                  // k is the row size and lda
                  DO I = N1, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, N1 - 1
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,n1+i)
                        WORK( I+N1 ) = WORK( I+N1 ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = S
                  END DO
                  // j=n1=k-1 is special
                  S = ABS( REAL( A( 0+J*LDA ) ) )
                  // A(k-1,k-1)
                  DO I = 1, K - 1
                     AA = ABS( A( I+J*LDA ) )
                     // A(k-1,i+n1)
                     WORK( I+N1 ) = WORK( I+N1 ) + AA
                     S = S + AA
                  END DO
                  WORK( J ) = WORK( J ) + S
                  DO J = K, N - 1
                     S = ZERO
                     DO I = 0, J - K - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(i,j-k)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     // i=j-k
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // A(j-k,j-k)
                     S = S + AA
                     WORK( J-K ) = WORK( J-K ) + S
                     I = I + 1
                     S = ABS( REAL( A( I+J*LDA ) ) )
                     // A(j,j)
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu=1 & uplo = 'L'
                  K = K + 1
                  // k=(n+1)/2 for n odd and ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 2
                     // process
                     S = ZERO
                     DO I = 0, J - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // i=j so process of A(j,j)
                     S = S + AA
                     WORK( J ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     S = AA
                     DO L = K + J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( K+J ) = WORK( K+J ) + S
                  END DO
                  // j=k-1 is special :process col A(k-1,0:k-1)
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO
                  // i=k-1
                  AA = ABS( REAL( A( I+J*LDA ) ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
                  // done with col j=k+1
                  DO J = K, N - 1
                     // process col j of A = A(j,0:k-1)
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            } else {
               // n is even & A is k=n/2 by n+1
               if ( ILU.EQ.0 ) {
                  // uplo = 'U'
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 1
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i+k)
                        WORK( I+K ) = WORK( I+K ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = S
                  END DO
                  // j=k
                  AA = ABS( REAL( A( 0+J*LDA ) ) )
                  // A(k,k)
                  S = AA
                  DO I = 1, K - 1
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,k+i)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  END DO
                  WORK( J ) = WORK( J ) + S
                  DO J = K + 1, N - 1
                     S = ZERO
                     DO I = 0, J - 2 - K
                        AA = ABS( A( I+J*LDA ) )
                        // A(i,j-k-1)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     // i=j-1-k
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // A(j-k-1,j-k-1)
                     S = S + AA
                     WORK( J-K-1 ) = WORK( J-K-1 ) + S
                     I = I + 1
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // A(j,j)
                     S = AA
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  // j=n
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
                     // A(i,k-1)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO
                  // i=k-1
                  AA = ABS( REAL( A( I+J*LDA ) ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = WORK( I ) + S
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu=1 & uplo = 'L'
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  // j=0 is special :process col A(k:n-1,k)
                  S = ABS( REAL( A( 0 ) ) )
                  // A(k,k)
                  DO I = 1, K - 1
                     AA = ABS( A( I ) )
                     // A(k+i,k)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  END DO
                  WORK( K ) = WORK( K ) + S
                  DO J = 1, K - 1
                     // process
                     S = ZERO
                     DO I = 0, J - 2
                        AA = ABS( A( I+J*LDA ) )
                        // A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     // i=j-1 so process of A(j-1,j-1)
                     S = S + AA
                     WORK( J-1 ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( REAL( A( I+J*LDA ) ) )
                     S = AA
                     DO L = K + J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( K+J ) = WORK( K+J ) + S
                  END DO
                  // j=k is special :process col A(k,0:k-1)
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO

                  // i=k-1
                  AA = ABS( REAL( A( I+J*LDA ) ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
                  // done with col j=k+1
                  DO J = K + 1, N

                     // process col j-1 of A = A(j-1,0:k-1)
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     WORK( J-1 ) = WORK( J-1 ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            }
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

        // Find normF(A).

         K = ( N+1 ) / 2
         SCALE = ZERO
         S = ONE
         if ( NOE.EQ.1 ) {
            // n is odd
            if ( IFM.EQ.1 ) {
               // A is normal & A is n by k
               if ( ILU.EQ.0 ) {
                  // A is upper
                  DO J = 0, K - 3
                     classq(K-J-2, A( K+J+1+J*LDA ), 1, SCALE, S );
                     // L at A(k,0)
                  END DO
                  DO J = 0, K - 1
                     classq(K+J-1, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = K - 1
                  // -> U(k,k) at A(k-1,0)
                  DO I = 0, K - 2
                     AA = REAL( A( L ) )
                     // U(k+i,k+i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // U(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
                  AA = REAL( A( L ) )
                  // U(n-1,n-1)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
               } else {
                  // ilu=1 & A is lower
                  DO J = 0, K - 1
                     classq(N-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // trap L at A(0,0)
                  END DO
                  DO J = 1, K - 2
                     classq(J, A( 0+( 1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  AA = REAL( A( 0 ) )
                  // L(0,0) at A(0,0)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
                  L = LDA
                  // -> L(k,k) at A(0,1)
                  DO I = 1, K - 1
                     AA = REAL( A( L ) )
                     // L(k-1+i,k-1+i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // L(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
               }
            } else {
               // A is xpose & A is k by n
               if ( ILU.EQ.0 ) {
                  // A**H is upper
                  DO J = 1, K - 2
                     classq(J, A( 0+( K+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k)
                  END DO
                  DO J = 0, K - 2
                     classq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     classq(K-J-1, A( J+1+( J+K-1 )*LDA ), 1, SCALE, S );
                     // L at A(0,k-1)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = 0 + K*LDA - LDA
                  // -> U(k-1,k-1) at A(0,k-1)
                  AA = REAL( A( L ) )
                  // U(k-1,k-1)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
                  L = L + LDA
                  // -> U(0,0) at A(0,k)
                  DO J = K, N - 1
                     AA = REAL( A( L ) )
                     // -> U(j-k,j-k)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // -> U(j,j)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
               } else {
                  // A**H is lower
                  DO J = 1, K - 1
                     classq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  END DO
                  DO J = K, N - 1
                     classq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,k)
                  END DO
                  DO J = 0, K - 3
                     classq(K-J-2, A( J+2+J*LDA ), 1, SCALE, S );
                     // L at A(1,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = 0
                  // -> L(0,0) at A(0,0)
                  DO I = 0, K - 2
                     AA = REAL( A( L ) )
                     // L(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // L(k+i,k+i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
                  // L-> k-1 + (k-1)*lda or L(k-1,k-1) at A(k-1,k-1)
                  AA = REAL( A( L ) )
                  // L(k-1,k-1) at A(k-1,k-1)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
               }
            }
         } else {
            // n is even
            if ( IFM.EQ.1 ) {
               // A is normal
               if ( ILU.EQ.0 ) {
                  // A is upper
                  DO J = 0, K - 2
                     classq(K-J-1, A( K+J+2+J*LDA ), 1, SCALE, S );
                  // L at A(k+1,0)
                  END DO
                  DO J = 0, K - 1
                     classq(K+J, A( 0+J*LDA ), 1, SCALE, S );
                  // trap U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = K
                  // -> U(k,k) at A(k,0)
                  DO I = 0, K - 1
                     AA = REAL( A( L ) )
                     // U(k+i,k+i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // U(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
               } else {
                  // ilu=1 & A is lower
                  DO J = 0, K - 1
                     classq(N-J-1, A( J+2+J*LDA ), 1, SCALE, S );
                     // trap L at A(1,0)
                  END DO
                  DO J = 1, K - 1
                     classq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = 0
                  // -> L(k,k) at A(0,0)
                  DO I = 0, K - 1
                     AA = REAL( A( L ) )
                     // L(k-1+i,k-1+i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // L(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
               }
            } else {
               // A is xpose
               if ( ILU.EQ.0 ) {
                  // A**H is upper
                  DO J = 1, K - 1
                     classq(J, A( 0+( K+1+J )*LDA ), 1, SCALE, S );
                  // U at A(0,k+1)
                  END DO
                  DO J = 0, K - 1
                     classq(K, A( 0+J*LDA ), 1, SCALE, S );
                  // k by k rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     classq(K-J-1, A( J+1+( J+K )*LDA ), 1, SCALE, S );
                  // L at A(0,k)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = 0 + K*LDA
                  // -> U(k,k) at A(0,k)
                  AA = REAL( A( L ) )
                  // U(k,k)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
                  L = L + LDA
                  // -> U(0,0) at A(0,k+1)
                  DO J = K + 1, N - 1
                     AA = REAL( A( L ) )
                     // -> U(j-k-1,j-k-1)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // -> U(j,j)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
                  // L=k-1+n*lda
                  // -> U(k-1,k-1) at A(k-1,n)
                  AA = REAL( A( L ) )
                  // U(k,k)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
               } else {
                  // A**H is lower
                  DO J = 1, K - 1
                     classq(J, A( 0+( J+1 )*LDA ), 1, SCALE, S );
                  // U at A(0,1)
                  END DO
                  DO J = K + 1, N
                     classq(K, A( 0+J*LDA ), 1, SCALE, S );
                  // k by k rect. at A(0,k+1)
                  END DO
                  DO J = 0, K - 2
                     classq(K-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                  // L at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  L = 0
                  // -> L(k,k) at A(0,0)
                  AA = REAL( A( L ) )
                  // L(k,k) at A(0,0)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
                  L = LDA
                  // -> L(0,0) at A(0,1)
                  DO I = 0, K - 2
                     AA = REAL( A( L ) )
                     // L(i,i)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     AA = REAL( A( L+1 ) )
                     // L(k+i+1,k+i+1)
                     if ( AA.NE.ZERO ) {
                        if ( SCALE.LT.AA ) {
                           S = ONE + S*( SCALE / AA )**2
                           SCALE = AA
                        } else {
                           S = S + ( AA / SCALE )**2
                        }
                     }
                     L = L + LDA + 1
                  END DO
                  // L-> k - 1 + k*lda or L(k-1,k-1) at A(k-1,k)
                  AA = REAL( A( L ) )
                  // L(k-1,k-1) at A(k-1,k)
                  if ( AA.NE.ZERO ) {
                     if ( SCALE.LT.AA ) {
                        S = ONE + S*( SCALE / AA )**2
                        SCALE = AA
                     } else {
                        S = S + ( AA / SCALE )**2
                     }
                  }
               }
            }
         }
         VALUE = SCALE*SQRT( S )
      }

      CLANHF = VALUE
      RETURN

      // End of CLANHF

      }
