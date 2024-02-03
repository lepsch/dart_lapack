      REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, TRANSR, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      REAL               A( 0: * ), WORK( 0: * )
      // ..

*  =====================================================================

      // ..
      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0, ZERO = 0.0 ;
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
      // EXTERNAL SLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         SLANSF = ZERO
         RETURN
      } else if ( N == 1 ) {
         SLANSF = ABS( A(0) )
         RETURN
      }

      // set noe = 1 if n is odd. if n is even set noe=0

      NOE = 1
      IF( MOD( N, 2 ) == 0 ) NOE = 0

      // set ifm = 0 when form='T or 't' and 1 otherwise

      IFM = 1
      IF( LSAME( TRANSR, 'T' ) ) IFM = 0

      // set ilu = 0 when uplo='U or 'u' and 1 otherwise

      ILU = 1
      IF( LSAME( UPLO, 'U' ) ) ILU = 0

      // set lda = (n+1)/2 when ifm = 0
      // set lda = n when ifm = 1 and noe = 1
      // set lda = n+1 when ifm = 1 and noe = 0

      if ( IFM == 1 ) {
         if ( NOE == 1 ) {
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
         if ( NOE == 1 ) {
            // n is odd
            if ( IFM == 1 ) {
            // A is n by k
               for (J = 0; J <= K - 1; J++) {
                  for (I = 0; I <= N - 1; I++) {
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            } else {
               // xpose case; A is k by n
               for (J = 0; J <= N - 1; J++) {
                  for (I = 0; I <= K - 1; I++) {
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            }
         } else {
            // n is even
            if ( IFM == 1 ) {
               // A is n+1 by k
               for (J = 0; J <= K - 1; J++) {
                  for (I = 0; I <= N; I++) {
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            } else {
               // xpose case; A is k by n+1
               for (J = 0; J <= N; J++) {
                  for (I = 0; I <= K - 1; I++) {
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            }
         }
      } else if ( ( LSAME( NORM, 'I' ) ) || ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         if ( IFM == 1 ) {
            K = N / 2
            if ( NOE == 1 ) {
               // n is odd
               if ( ILU == 0 ) {
                  for (I = 0; I <= K - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  for (J = 0; J <= K; J++) {
                     S = ZERO
                     for (I = 0; I <= K + J - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     }
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     if (I == K+K) GO TO 10;
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     for (L = J + 1; L <= K - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  } // 10
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               } else {
                  // ilu = 1
                  K = K + 1
                  // k=(n+1)/2 for n odd and ilu=1
                  for (I = K; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  DO J = K - 1, 0, -1
                     S = ZERO
                     for (I = 0; I <= J - 2; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     }
                     if ( J > 0 ) {
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
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            } else {
               // n is even
               if ( ILU == 0 ) {
                  for (I = 0; I <= K - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  for (J = 0; J <= K - 1; J++) {
                     S = ZERO
                     for (I = 0; I <= K + J - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     }
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
                     // -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     for (L = J + 1; L <= K - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               } else {
                  // ilu = 1
                  for (I = K; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  DO J = K - 1, 0, -1
                     S = ZERO
                     for (I = 0; I <= J - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     }
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
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            }
         } else {
            // ifm=0
            K = N / 2
            if ( NOE == 1 ) {
               // n is odd
               if ( ILU == 0 ) {
                  N1 = K
                  // n/2
                  K = K + 1
                  // k is the row size and lda
                  for (I = N1; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  for (J = 0; J <= N1 - 1; J++) {
                     S = ZERO
                     for (I = 0; I <= K - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,n1+i)
                        WORK( I+N1 ) = WORK( I+N1 ) + AA
                        S = S + AA
                     }
                     WORK( J ) = S
                  }
                  // j=n1=k-1 is special
                  S = ABS( A( 0+J*LDA ) )
                  // A(k-1,k-1)
                  for (I = 1; I <= K - 1; I++) {
                     AA = ABS( A( I+J*LDA ) )
                     // A(k-1,i+n1)
                     WORK( I+N1 ) = WORK( I+N1 ) + AA
                     S = S + AA
                  }
                  WORK( J ) = WORK( J ) + S
                  for (J = K; J <= N - 1; J++) {
                     S = ZERO
                     for (I = 0; I <= J - K - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(i,j-k)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                     // i=j-k
                     AA = ABS( A( I+J*LDA ) )
                     // A(j-k,j-k)
                     S = S + AA
                     WORK( J-K ) = WORK( J-K ) + S
                     I = I + 1
                     S = ABS( A( I+J*LDA ) )
                     // A(j,j)
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               } else {
                  // ilu=1
                  K = K + 1
                  // k=(n+1)/2 for n odd and ilu=1
                  for (I = K; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  for (J = 0; J <= K - 2; J++) {
                     // process
                     S = ZERO
                     for (I = 0; I <= J - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                     AA = ABS( A( I+J*LDA ) )
                     // i=j so process of A(j,j)
                     S = S + AA
                     WORK( J ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
                     S = AA
                     for (L = K + J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( K+J ) = WORK( K+J ) + S
                  }
                  // j=k-1 is special :process col A(k-1,0:k-1)
                  S = ZERO
                  for (I = 0; I <= K - 2; I++) {
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  }
                  // i=k-1
                  AA = ABS( A( I+J*LDA ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
                  // done with col j=k+1
                  for (J = K; J <= N - 1; J++) {
                     // process col j of A = A(j,0:k-1)
                     S = ZERO
                     for (I = 0; I <= K - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            } else {
               // n is even
               if ( ILU == 0 ) {
                  for (I = K; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  for (J = 0; J <= K - 1; J++) {
                     S = ZERO
                     for (I = 0; I <= K - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,i+k)
                        WORK( I+K ) = WORK( I+K ) + AA
                        S = S + AA
                     }
                     WORK( J ) = S
                  }
                  // j=k
                  AA = ABS( A( 0+J*LDA ) )
                  // A(k,k)
                  S = AA
                  for (I = 1; I <= K - 1; I++) {
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,k+i)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  }
                  WORK( J ) = WORK( J ) + S
                  for (J = K + 1; J <= N - 1; J++) {
                     S = ZERO
                     for (I = 0; I <= J - 2 - K; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(i,j-k-1)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                      // i=j-1-k
                     AA = ABS( A( I+J*LDA ) )
                     // A(j-k-1,j-k-1)
                     S = S + AA
                     WORK( J-K-1 ) = WORK( J-K-1 ) + S
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
                     // A(j,j)
                     S = AA
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     }
                     WORK( J ) = WORK( J ) + S
                  }
                  // j=n
                  S = ZERO
                  for (I = 0; I <= K - 2; I++) {
                     AA = ABS( A( I+J*LDA ) )
                     // A(i,k-1)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  }
                  // i=k-1
                  AA = ABS( A( I+J*LDA ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = WORK( I ) + S
                  VALUE = WORK ( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               } else {
                  // ilu=1
                  for (I = K; I <= N - 1; I++) {
                     WORK( I ) = ZERO
                  }
                  // j=0 is special :process col A(k:n-1,k)
                  S = ABS( A( 0 ) )
                  // A(k,k)
                  for (I = 1; I <= K - 1; I++) {
                     AA = ABS( A( I ) )
                     // A(k+i,k)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  }
                  WORK( K ) = WORK( K ) + S
                  for (J = 1; J <= K - 1; J++) {
                     // process
                     S = ZERO
                     for (I = 0; I <= J - 2; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                     AA = ABS( A( I+J*LDA ) )
                     // i=j-1 so process of A(j-1,j-1)
                     S = S + AA
                     WORK( J-1 ) = S
                     // is initialised here
                     I = I + 1
                     // i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
                     S = AA
                     for (L = K + J + 1; L <= N - 1; L++) {
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
                        // A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     }
                     WORK( K+J ) = WORK( K+J ) + S
                  }
                  // j=k is special :process col A(k,0:k-1)
                  S = ZERO
                  for (I = 0; I <= K - 2; I++) {
                     AA = ABS( A( I+J*LDA ) )
                     // A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  }
                  // i=k-1
                  AA = ABS( A( I+J*LDA ) )
                  // A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
                  // done with col j=k+1
                  for (J = K + 1; J <= N; J++) {
                     // process col j-1 of A = A(j-1,0:k-1)
                     S = ZERO
                     for (I = 0; I <= K - 1; I++) {
                        AA = ABS( A( I+J*LDA ) )
                        // A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     }
                     WORK( J-1 ) = WORK( J-1 ) + S
                  }
                  VALUE = WORK( 0 )
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I )
                     IF( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP
                  }
               }
            }
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

        // Find normF(A).

         K = ( N+1 ) / 2
         SCALE = ZERO
         S = ONE
         if ( NOE == 1 ) {
            // n is odd
            if ( IFM == 1 ) {
               // A is normal
               if ( ILU == 0 ) {
                  // A is upper
                  for (J = 0; J <= K - 3; J++) {
                     slassq(K-J-2, A( K+J+1+J*LDA ), 1, SCALE, S );
                     // L at A(k,0)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     slassq(K+J-1, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K-1, A( K ), LDA+1, SCALE, S );
                  // tri L at A(k,0)
                  slassq(K, A( K-1 ), LDA+1, SCALE, S );
                  // tri U at A(k-1,0)
               } else {
                  // ilu=1 & A is lower
                  for (J = 0; J <= K - 1; J++) {
                     slassq(N-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // trap L at A(0,0)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     slassq(J, A( 0+( 1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri L at A(0,0)
                  slassq(K-1, A( 0+LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,1)
               }
            } else {
               // A is xpose
               if ( ILU == 0 ) {
                  // A**T is upper
                  for (J = 1; J <= K - 2; J++) {
                     slassq(J, A( 0+( K+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     slassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,0)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     slassq(K-J-1, A( J+1+( J+K-1 )*LDA ), 1, SCALE, S );
                     // L at A(0,k-1)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K-1, A( 0+K*LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,k)
                  slassq(K, A( 0+( K-1 )*LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,k-1)
               } else {
                  // A**T is lower
                  for (J = 1; J <= K - 1; J++) {
                     slassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  }
                  for (J = K; J <= N - 1; J++) {
                     slassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,k)
                  }
                  for (J = 0; J <= K - 3; J++) {
                     slassq(K-J-2, A( J+2+J*LDA ), 1, SCALE, S );
                     // L at A(1,0)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
                  slassq(K-1, A( 1 ), LDA+1, SCALE, S );
                  // tri L at A(1,0)
               }
            }
         } else {
            // n is even
            if ( IFM == 1 ) {
               // A is normal
               if ( ILU == 0 ) {
                  // A is upper
                  for (J = 0; J <= K - 2; J++) {
                     slassq(K-J-1, A( K+J+2+J*LDA ), 1, SCALE, S );
                     // L at A(k+1,0)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     slassq(K+J, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( K+1 ), LDA+1, SCALE, S );
                  // tri L at A(k+1,0)
                  slassq(K, A( K ), LDA+1, SCALE, S );
                  // tri U at A(k,0)
               } else {
                  // ilu=1 & A is lower
                  for (J = 0; J <= K - 1; J++) {
                     slassq(N-J-1, A( J+2+J*LDA ), 1, SCALE, S );
                     // trap L at A(1,0)
                  }
                  for (J = 1; J <= K - 1; J++) {
                     slassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( 1 ), LDA+1, SCALE, S );
                  // tri L at A(1,0)
                  slassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
               }
            } else {
               // A is xpose
               if ( ILU == 0 ) {
                  // A**T is upper
                  for (J = 1; J <= K - 1; J++) {
                     slassq(J, A( 0+( K+1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k+1)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     slassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k rect. at A(0,0)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     slassq(K-J-1, A( J+1+( J+K )*LDA ), 1, SCALE, S );
                     // L at A(0,k)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( 0+( K+1 )*LDA ), LDA+1, SCALE, S );
                  // tri U at A(0,k+1)
                  slassq(K, A( 0+K*LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,k)
               } else {
                  // A**T is lower
                  for (J = 1; J <= K - 1; J++) {
                     slassq(J, A( 0+( J+1 )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  }
                  for (J = K + 1; J <= N; J++) {
                     slassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k rect. at A(0,k+1)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     slassq(K-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // L at A(0,0)
                  }
                  S = S + S
                  // double s for the off diagonal elements
                  slassq(K, A( LDA ), LDA+1, SCALE, S );
                  // tri L at A(0,1)
                  slassq(K, A( 0 ), LDA+1, SCALE, S );
                  // tri U at A(0,0)
               }
            }
         }
         VALUE = SCALE*SQRT( S )
      }

      SLANSF = VALUE
      RETURN

      // End of SLANSF

      }
