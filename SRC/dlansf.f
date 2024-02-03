      double           FUNCTION DLANSF( NORM, TRANSR, UPLO, N, A, WORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, TRANSR, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      double             A( 0: * ), WORK( 0: * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, IFM, ILU, NOE, N1, K, L, LDA;
      double             SCALE, S, VALUE, AA, TEMP;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         DLANSF = ZERO
         RETURN
      } else if ( N.EQ.1 ) {
         DLANSF = ABS( A(0) )
         RETURN
      }

      // set noe = 1 if n is odd. if n is even set noe=0

      NOE = 1
      IF( MOD( N, 2 ).EQ.0 ) NOE = 0

      // set ifm = 0 when form='T or 't' and 1 otherwise

      IFM = 1
      IF( LSAME( TRANSR, 'T' ) ) IFM = 0

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
            // n is odd
            if ( IFM.EQ.1 ) {
            // A is n by k
               DO J = 0, K - 1
                  DO I = 0, N - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               END DO
            } else {
               // xpose case; A is k by n
               DO J = 0, N - 1
                  DO I = 0, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               END DO
            }
         } else {
            // n is even
            if ( IFM.EQ.1 ) {
               // A is n+1 by k
               DO J = 0, K - 1
                  DO I = 0, N
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               END DO
            } else {
               // xpose case; A is k by n+1
               DO J = 0, N
                  DO I = 0, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               END DO
            }
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         if ( IFM.EQ.1 ) {
            K = N / 2
            if ( NOE.EQ.1 ) {
               // n is odd
               if ( ILU.EQ.0 ) {
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
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     IF( I.EQ.K+K ) GO TO 10
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu = 1
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
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(j+k,j+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + S
                        // i=j
                        I = I + 1
                     }
                     AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            } else {
               // n is even
               if ( ILU.EQ.0 ) {
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
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu = 1
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
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j+k,j+k)
                     S = S + AA
                     WORK( I+K ) = WORK( I+K ) + S
                     // i=j
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            }
         } else {
            // ifm=0
            K = N / 2
            if ( NOE.EQ.1 ) {
               // n is odd
               if ( ILU.EQ.0 ) {
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
                  S = ABS( A( 0+J*LDA ) )
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
                     AA = ABS( A( I+J*LDA ) )
                     // A(j-k,j-k)
                     S = S + AA
                     WORK( J-K ) = WORK( J-K ) + S
                     I = I + 1
                     S = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu=1
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
                     AA = ABS( A( I+J*LDA ) )
                     // i=j so process of A(j,j)
                     S = S + AA
                     WORK( J ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
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
                  AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               }
            } else {
               // n is even
               if ( ILU.EQ.0 ) {
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
                  AA = ABS( A( 0+J*LDA ) )
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
                     AA = ABS( A( I+J*LDA ) )
                     // A(j-k-1,j-k-1)
                     S = S + AA
                     WORK( J-K-1 ) = WORK( J-K-1 ) + S
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
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
                  AA = ABS( A( I+J*LDA ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = WORK( I ) + S
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
                  END DO
               } else {
                  // ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  // j=0 is special :process col A(k:n-1,k)
                  S = ABS( A( 0 ) )
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
                     AA = ABS( A( I+J*LDA ) )
                     // i=j-1 so process of A(j-1,j-1)
                     S = S + AA
                     WORK( J-1 ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
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
                  AA = ABS( A( I+J*LDA ) )
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
                     IF( VALUE .LT. TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
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
               // A is normal
               if ( ILU.EQ.0 ) {
                  // A is upper
                  DO J = 0, K - 3
                     dlassq(K-J-2, A( K+J+1+J*LDA ), 1, SCALE, S );
                     // L at A(k,0)
                  END DO
                  DO J = 0, K - 1
                     dlassq(K+J-1, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K-1, A( K ), LDA+1, SCALE, S );
                  // tri L at A(k,0)
                  dlassq(K, A( K-1 ), LDA+1, SCALE, S );
                  // tri U at A(k-1,0)
               } else {
                  // ilu=1 & A is lower
                  DO J = 0, K - 1
                     dlassq(N-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // trap L at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     dlassq(J, A( 0+( 1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri L at A(0,0)
                  dlassq(K-1, A( 0+LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,1)
               }
            } else {
               // A is xpose
               if ( ILU.EQ.0 ) {
                  // A**T is upper
                  DO J = 1, K - 2
                     dlassq(J, A( 0+( K+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k)
                  END DO
                  DO J = 0, K - 2
                     dlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     dlassq(K-J-1, A( J+1+( J+K-1 )*LDA ), 1, SCALE, S );
                     // L at A(0,k-1)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K-1, A( 0+K*LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,k)
                  dlassq(K, A( 0+( K-1 )*LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,k-1)
               } else {
                  // A**T is lower
                  DO J = 1, K - 1
                     dlassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  END DO
                  DO J = K, N - 1
                     dlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,k)
                  END DO
                  DO J = 0, K - 3
                     dlassq(K-J-2, A( J+2+J*LDA ), 1, SCALE, S );
                     // L at A(1,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
                  dlassq(K-1, A( 1 ), LDA+1, SCALE, S );
                  // tri L at A(1,0)
               }
            }
         } else {
            // n is even
            if ( IFM.EQ.1 ) {
               // A is normal
               if ( ILU.EQ.0 ) {
                  // A is upper
                  DO J = 0, K - 2
                     dlassq(K-J-1, A( K+J+2+J*LDA ), 1, SCALE, S );
                     // L at A(k+1,0)
                  END DO
                  DO J = 0, K - 1
                     dlassq(K+J, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( K+1 ), LDA+1, SCALE, S );
                  // tri L at A(k+1,0)
                  dlassq(K, A( K ), LDA+1, SCALE, S );
                  // tri U at A(k,0)
               } else {
                  // ilu=1 & A is lower
                  DO J = 0, K - 1
                     dlassq(N-J-1, A( J+2+J*LDA ), 1, SCALE, S );
                     // trap L at A(1,0)
                  END DO
                  DO J = 1, K - 1
                     dlassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( 1 ), LDA+1, SCALE, S );
                  // tri L at A(1,0)
                  dlassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
               }
            } else {
               // A is xpose
               if ( ILU.EQ.0 ) {
                  // A**T is upper
                  DO J = 1, K - 1
                     dlassq(J, A( 0+( K+1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k+1)
                  END DO
                  DO J = 0, K - 1
                     dlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     dlassq(K-J-1, A( J+1+( J+K )*LDA ), 1, SCALE, S );
                     // L at A(0,k)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( 0+( K+1 )*LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,k+1)
                  dlassq(K, A( 0+K*LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,k)
               } else {
                  // A**T is lower
                  DO J = 1, K - 1
                     dlassq(J, A( 0+( J+1 )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  END DO
                  DO J = K + 1, N
                     dlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k rect. at A(0,k+1)
                  END DO
                  DO J = 0, K - 2
                     dlassq(K-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // L at A(0,0)
                  END DO
                  S = S + S
                  // double s for the off diagonal elements
                  dlassq(K, A( LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,1)
                  dlassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
               }
            }
         }
         VALUE = SCALE*SQRT( S )
      }

      DLANSF = VALUE
      RETURN

      // End of DLANSF

      }
