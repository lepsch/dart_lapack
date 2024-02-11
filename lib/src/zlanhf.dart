      double zlanhf(final int NORM, final int TRANSR, final int UPLO, final int N, final int A, final Array<double> WORK,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM, TRANSR, UPLO;
      int                N;
      double             WORK( 0: * );
      Complex         A( 0: * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J, IFM, ILU, NOE, N1, K, L, LDA;
      double             SCALE, S, VALUE, AA, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      // EXTERNAL lsame, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, SQRT

      if ( N == 0 ) {
         ZLANHF = ZERO;
         return;
      } else if ( N == 1 ) {
         ZLANHF = ABS((A(0))).toDouble();
         return;
      }

      // set noe = 1 if n is odd. if n is even set noe=0

      NOE = 1;
      if( (N % 2) == 0 ) NOE = 0;

      // set ifm = 0 when form='C' or 'c' and 1 otherwise

      IFM = 1;
      if( lsame( TRANSR, 'C' ) ) IFM = 0;

      // set ilu = 0 when uplo='U or 'u' and 1 otherwise

      ILU = 1;
      if( lsame( UPLO, 'U' ) ) ILU = 0;

      // set lda = (n+1)/2 when ifm = 0
      // set lda = n when ifm = 1 and noe = 1
      // set lda = n+1 when ifm = 1 and noe = 0

      if ( IFM == 1 ) {
         if ( NOE == 1 ) {
            LDA = N;
         } else {
            // noe=0
            LDA = N + 1;
         }
      } else {
         // ifm=0
         LDA = ( N+1 ) / 2;
      }

      if ( lsame( NORM, 'M' ) ) {

        // Find max(abs(A(i,j))).

         K = ( N+1 ) / 2;
         VALUE = ZERO;
         if ( NOE == 1 ) {
            // n is odd & n = k + k - 1
            if ( IFM == 1 ) {
               // A is n by k
               if ( ILU == 1 ) {
                  // uplo ='L'
                  J = 0;
                  // -> L(0,0)
                  TEMP = ABS( (A( J+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (I = 1; I <= N - 1; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  for (J = 1; J <= K - 1; J++) {
                     for (I = 0; I <= J - 2; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J - 1;
                     // L(k+j,k+j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J;
                     // -> L(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J + 1; I <= N - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
               } else {
                  // uplo = 'U'
                  for (J = 0; J <= K - 2; J++) {
                     for (I = 0; I <= K + J - 2; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = K + J - 1;
                     // -> U(i,i)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = I + 1;
                     // =k+j; i -> U(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = K + J + 1; I <= N - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  for (I = 0; I <= N - 2; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     // j=k-1
                  }
                  // i=n-1 -> U(n-1,n-1)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
               }
            } else {
               // xpose case; A is k by n
               if ( ILU == 1 ) {
                  // uplo ='L'
                  for (J = 0; J <= K - 2; J++) {
                     for (I = 0; I <= J - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J;
                     // L(i,i)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J + 1;
                     // L(j+k,j+k)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J + 2; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  J = K - 1;
                  for (I = 0; I <= K - 2; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  I = K - 1;
                  // -> L(i,i) is at A(i,j)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (J = K; J <= N - 1; J++) {
                     for (I = 0; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
               } else {
                  // uplo = 'U'
                  for (J = 0; J <= K - 2; J++) {
                     for (I = 0; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  J = K - 1;
                  // -> U(j,j) is at A(0,j)
                  TEMP = ABS( (A( 0+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (I = 1; I <= K - 1; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  for (J = K; J <= N - 1; J++) {
                     for (I = 0; I <= J - K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J - K;
                     // -> U(i,i) at A(i,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J - K + 1;
                     // U(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J - K + 2; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
               }
            }
         } else {
            // n is even & k = n/2
            if ( IFM == 1 ) {
               // A is n+1 by k
               if ( ILU == 1 ) {
                  // uplo ='L'
                  J = 0;
                  // -> L(k,k) & j=1 -> L(0,0)
                  TEMP = ABS( (A( J+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  TEMP = ABS( (A( J+1+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (I = 2; I <= N; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  for (J = 1; J <= K - 1; J++) {
                     for (I = 0; I <= J - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J;
                     // L(k+j,k+j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J + 1;
                     // -> L(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J + 2; I <= N; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
               } else {
                  // uplo = 'U'
                  for (J = 0; J <= K - 2; J++) {
                     for (I = 0; I <= K + J - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = K + J;
                     // -> U(i,i)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = I + 1;
                     // =k+j+1; i -> U(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = K + J + 2; I <= N; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  for (I = 0; I <= N - 2; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     // j=k-1
                  }
                  // i=n-1 -> U(n-1,n-1)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  I = N;
                  // -> U(k-1,k-1)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
               }
            } else {
               // xpose case; A is k by n+1
               if ( ILU == 1 ) {
                  // uplo ='L'
                  J = 0;
                  // -> L(k,k) at A(0,0)
                  TEMP = ABS( (A( J+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (I = 1; I <= K - 1; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  for (J = 1; J <= K - 1; J++) {
                     for (I = 0; I <= J - 2; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J - 1;
                     // L(i,i)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J;
                     // L(j+k,j+k)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J + 1; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  J = K;
                  for (I = 0; I <= K - 2; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  I = K - 1;
                  // -> L(i,i) is at A(i,j)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (J = K + 1; J <= N; J++) {
                     for (I = 0; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
               } else {
                  // uplo = 'U'
                  for (J = 0; J <= K - 1; J++) {
                     for (I = 0; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  J = K;
                  // -> U(j,j) is at A(0,j)
                  TEMP = ABS( (A( 0+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  for (I = 1; I <= K - 1; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  for (J = K + 1; J <= N - 1; J++) {
                     for (I = 0; I <= J - K - 2; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                     I = J - K - 1;
                     // -> U(i,i) at A(i,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     I = J - K;
                     // U(j,j)
                     TEMP = ABS( (A( I+J*LDA )).toDouble() );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     for (I = J - K + 1; I <= K - 1; I++) {
                        TEMP = ( A( I+J*LDA ) ).abs();
                        if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                     }
                  }
                  J = N;
                  for (I = 0; I <= K - 2; I++) {
                     TEMP = ( A( I+J*LDA ) ).abs();
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
                  I = K - 1;
                  // U(k,k) at A(i,j)
                  TEMP = ABS( (A( I+J*LDA )).toDouble() );
                  if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
               }
            }
         }
      } else if ( ( lsame( NORM, 'I' ) ) || ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

        // Find normI(A) ( = norm1(A), since A is Hermitian).

         if ( IFM == 1 ) {
            // A is 'N'
            K = N / 2;
            if ( NOE == 1 ) {
               // n is odd & A is n by (n+1)/2
               if ( ILU == 0 ) {
                  // uplo = 'U'
                  for (I = 0; I <= K - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = 0; J <= K; J++) {
                     S = ZERO;
                     for (I = 0; I <= K + J - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(i,j+k)
                        S = S + AA;
                        WORK[I] = WORK( I ) + AA;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j+k,j+k)
                     WORK[J+K] = S + AA;
                     if (I == K+K) GO TO 10;
                     I = I + 1;
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j,j)
                     WORK[J] = WORK( J ) + AA;
                     S = ZERO;
                     for (L = J + 1; L <= K - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(l,j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  } // 10
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               } else {
                  // ilu = 1 & uplo = 'L'
                  K = K + 1;
                  // k=(n+1)/2 for n odd and ilu=1
                  for (I = K; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = K - 1; J >= 0; J--) {
                     S = ZERO;
                     for (I = 0; I <= J - 2; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(j+k,i+k)
                        S = S + AA;
                        WORK[I+K] = WORK( I+K ) + AA;
                     }
                     if ( J > 0 ) {
                        AA = ABS( (A( I+J*LDA )).toDouble() );
                        // -> A(j+k,j+k)
                        S = S + AA;
                        WORK[I+K] = WORK( I+K ) + S;
                        // i=j
                        I = I + 1;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j,j)
                     WORK[J] = AA;
                     S = ZERO;
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(l,j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               }
            } else {
               // n is even & A is n+1 by k = n/2
               if ( ILU == 0 ) {
                  // uplo = 'U'
                  for (I = 0; I <= K - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = 0; J <= K - 1; J++) {
                     S = ZERO;
                     for (I = 0; I <= K + J - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(i,j+k)
                        S = S + AA;
                        WORK[I] = WORK( I ) + AA;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j+k,j+k)
                     WORK[J+K] = S + AA;
                     I = I + 1;
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j,j)
                     WORK[J] = WORK( J ) + AA;
                     S = ZERO;
                     for (L = J + 1; L <= K - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(l,j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               } else {
                  // ilu = 1 & uplo = 'L'
                  for (I = K; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = K - 1; J >= 0; J--) {
                     S = ZERO;
                     for (I = 0; I <= J - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(j+k,i+k)
                        S = S + AA;
                        WORK[I+K] = WORK( I+K ) + AA;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j+k,j+k)
                     S = S + AA;
                     WORK[I+K] = WORK( I+K ) + S;
                     // i=j
                     I = I + 1;
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // -> A(j,j)
                     WORK[J] = AA;
                     S = ZERO;
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // -> A(l,j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               }
            }
         } else {
            // ifm=0
            K = N / 2;
            if ( NOE == 1 ) {
               // n is odd & A is (n+1)/2 by n
               if ( ILU == 0 ) {
                  // uplo = 'U'
                  N1 = K;
                  // n/2
                  K = K + 1;
                  // k is the row size and lda
                  for (I = N1; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = 0; J <= N1 - 1; J++) {
                     S = ZERO;
                     for (I = 0; I <= K - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,n1+i)
                        WORK[I+N1] = WORK( I+N1 ) + AA;
                        S = S + AA;
                     }
                     WORK[J] = S;
                  }
                  // j=n1=k-1 is special
                  S = ABS( (A( 0+J*LDA )).toDouble() );
                  // A(k-1,k-1)
                  for (I = 1; I <= K - 1; I++) {
                     AA = ( A( I+J*LDA ) ).abs();
                     // A(k-1,i+n1)
                     WORK[I+N1] = WORK( I+N1 ) + AA;
                     S = S + AA;
                  }
                  WORK[J] = WORK( J ) + S;
                  for (J = K; J <= N - 1; J++) {
                     S = ZERO;
                     for (I = 0; I <= J - K - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(i,j-k)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     // i=j-k
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // A(j-k,j-k)
                     S = S + AA;
                     WORK[J-K] = WORK( J-K ) + S;
                     I = I + 1;
                     S = ABS( (A( I+J*LDA )).toDouble() );
                     // A(j,j)
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,l)
                        WORK[L] = WORK( L ) + AA;
                        S = S + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               } else {
                  // ilu=1 & uplo = 'L'
                  K = K + 1;
                  // k=(n+1)/2 for n odd and ilu=1
                  for (I = K; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = 0; J <= K - 2; J++) {
                     // process
                     S = ZERO;
                     for (I = 0; I <= J - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,i)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // i=j so process of A(j,j)
                     S = S + AA;
                     WORK[J] = S;
                     // is initialised here
                     I = I + 1;
                     // i=j process A(j+k,j+k)
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     S = AA;
                     for (L = K + J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(l,k+j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[K+J] = WORK( K+J ) + S;
                  }
                  // j=k-1 is special :process col A(k-1,0:k-1)
                  S = ZERO;
                  for (I = 0; I <= K - 2; I++) {
                     AA = ( A( I+J*LDA ) ).abs();
                     // A(k,i)
                     WORK[I] = WORK( I ) + AA;
                     S = S + AA;
                  }
                  // i=k-1
                  AA = ABS( (A( I+J*LDA )).toDouble() );
                  // A(k-1,k-1)
                  S = S + AA;
                  WORK[I] = S;
                  // done with col j=k+1
                  for (J = K; J <= N - 1; J++) {
                     // process col j of A = A(j,0:k-1)
                     S = ZERO;
                     for (I = 0; I <= K - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,i)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               }
            } else {
               // n is even & A is k=n/2 by n+1
               if ( ILU == 0 ) {
                  // uplo = 'U'
                  for (I = K; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  for (J = 0; J <= K - 1; J++) {
                     S = ZERO;
                     for (I = 0; I <= K - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,i+k)
                        WORK[I+K] = WORK( I+K ) + AA;
                        S = S + AA;
                     }
                     WORK[J] = S;
                  }
                  // j=k
                  AA = ABS( (A( 0+J*LDA )).toDouble() );
                  // A(k,k)
                  S = AA;
                  for (I = 1; I <= K - 1; I++) {
                     AA = ( A( I+J*LDA ) ).abs();
                     // A(k,k+i)
                     WORK[I+K] = WORK( I+K ) + AA;
                     S = S + AA;
                  }
                  WORK[J] = WORK( J ) + S;
                  for (J = K + 1; J <= N - 1; J++) {
                     S = ZERO;
                     for (I = 0; I <= J - 2 - K; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(i,j-k-1)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     // i=j-1-k
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // A(j-k-1,j-k-1)
                     S = S + AA;
                     WORK[J-K-1] = WORK( J-K-1 ) + S;
                     I = I + 1;
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // A(j,j)
                     S = AA;
                     for (L = J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j,l)
                        WORK[L] = WORK( L ) + AA;
                        S = S + AA;
                     }
                     WORK[J] = WORK( J ) + S;
                  }
                  // j=n
                  S = ZERO;
                  for (I = 0; I <= K - 2; I++) {
                     AA = ( A( I+J*LDA ) ).abs();
                     // A(i,k-1)
                     WORK[I] = WORK( I ) + AA;
                     S = S + AA;
                  }
                  // i=k-1
                  AA = ABS( (A( I+J*LDA )).toDouble() );
                  // A(k-1,k-1)
                  S = S + AA;
                  WORK[I] = WORK( I ) + S;
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               } else {
                  // ilu=1 & uplo = 'L'
                  for (I = K; I <= N - 1; I++) {
                     WORK[I] = ZERO;
                  }
                  // j=0 is special :process col A(k:n-1,k)
                  S = ABS( (A( 0 )).toDouble() );
                  // A(k,k)
                  for (I = 1; I <= K - 1; I++) {
                     AA = ( A( I ) ).abs();
                     // A(k+i,k)
                     WORK[I+K] = WORK( I+K ) + AA;
                     S = S + AA;
                  }
                  WORK[K] = WORK( K ) + S;
                  for (J = 1; J <= K - 1; J++) {
                     // process
                     S = ZERO;
                     for (I = 0; I <= J - 2; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j-1,i)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     // i=j-1 so process of A(j-1,j-1)
                     S = S + AA;
                     WORK[J-1] = S;
                     // is initialised here
                     I = I + 1;
                     // i=j process A(j+k,j+k)
                     AA = ABS( (A( I+J*LDA )).toDouble() );
                     S = AA;
                     for (L = K + J + 1; L <= N - 1; L++) {
                        I = I + 1;
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(l,k+j)
                        S = S + AA;
                        WORK[L] = WORK( L ) + AA;
                     }
                     WORK[K+J] = WORK( K+J ) + S;
                  }
                  // j=k is special :process col A(k,0:k-1)
                  S = ZERO;
                  for (I = 0; I <= K - 2; I++) {
                     AA = ( A( I+J*LDA ) ).abs();
                     // A(k,i)
                     WORK[I] = WORK( I ) + AA;
                     S = S + AA;
                  }

                  // i=k-1
                  AA = ABS( (A( I+J*LDA )).toDouble() );
                  // A(k-1,k-1)
                  S = S + AA;
                  WORK[I] = S;
                  // done with col j=k+1
                  for (J = K + 1; J <= N; J++) {

                     // process col j-1 of A = A(j-1,0:k-1)
                     S = ZERO;
                     for (I = 0; I <= K - 1; I++) {
                        AA = ( A( I+J*LDA ) ).abs();
                        // A(j-1,i)
                        WORK[I] = WORK( I ) + AA;
                        S = S + AA;
                     }
                     WORK[J-1] = WORK( J-1 ) + S;
                  }
                  VALUE = WORK( 0 );
                  for (I = 1; I <= N-1; I++) {
                     TEMP = WORK( I );
                     if( VALUE < TEMP || disnan( TEMP ) ) VALUE = TEMP;
                  }
               }
            }
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

        // Find normF(A).

         K = ( N+1 ) / 2;
         SCALE = ZERO;
         S = ONE;
         if ( NOE == 1 ) {
            // n is odd
            if ( IFM == 1 ) {
               // A is normal & A is n by k
               if ( ILU == 0 ) {
                  // A is upper
                  for (J = 0; J <= K - 3; J++) {
                     zlassq(K-J-2, A( K+J+1+J*LDA ), 1, SCALE, S );
                     // L at A(k,0)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     zlassq(K+J-1, A( 0+J*LDA ), 1, SCALE, S );
                     // trap U at A(0,0)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = K - 1;
                  // -> U(k,k) at A(k-1,0)
                  for (I = 0; I <= K - 2; I++) {
                     AA = (A( L )).toDouble();
                     // U(k+i,k+i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // U(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
                  AA = (A( L )).toDouble();
                  // U(n-1,n-1)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
               } else {
                  // ilu=1 & A is lower
                  for (J = 0; J <= K - 1; J++) {
                     zlassq(N-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                     // trap L at A(0,0)
                  }
                  for (J = 1; J <= K - 2; J++) {
                     zlassq(J, A( 0+( 1+J )*LDA ), 1, SCALE, S );
                     // U at A(0,1)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  AA = (A( 0 )).toDouble();
                  // L(0,0) at A(0,0)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
                  L = LDA;
                  // -> L(k,k) at A(0,1)
                  for (I = 1; I <= K - 1; I++) {
                     AA = (A( L )).toDouble();
                     // L(k-1+i,k-1+i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // L(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
               }
            } else {
               // A is xpose & A is k by n
               if ( ILU == 0 ) {
                  // A**H is upper
                  for (J = 1; J <= K - 2; J++) {
                     zlassq(J, A( 0+( K+J )*LDA ), 1, SCALE, S );
                     // U at A(0,k)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     zlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,0)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     zlassq(K-J-1, A( J+1+( J+K-1 )*LDA ), 1, SCALE, S );
                     // L at A(0,k-1)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = 0 + K*LDA - LDA;
                  // -> U(k-1,k-1) at A(0,k-1)
                  AA = (A( L )).toDouble();
                  // U(k-1,k-1)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
                  L = L + LDA;
                  // -> U(0,0) at A(0,k)
                  for (J = K; J <= N - 1; J++) {
                     AA = (A( L )).toDouble();
                     // -> U(j-k,j-k)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // -> U(j,j)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
               } else {
                  // A**H is lower
                  for (J = 1; J <= K - 1; J++) {
                     zlassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  }
                  for (J = K; J <= N - 1; J++) {
                     zlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                     // k by k-1 rect. at A(0,k)
                  }
                  for (J = 0; J <= K - 3; J++) {
                     zlassq(K-J-2, A( J+2+J*LDA ), 1, SCALE, S );
                     // L at A(1,0)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = 0;
                  // -> L(0,0) at A(0,0)
                  for (I = 0; I <= K - 2; I++) {
                     AA = (A( L )).toDouble();
                     // L(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // L(k+i,k+i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
                  // L-> k-1 + (k-1)*lda or L(k-1,k-1) at A(k-1,k-1)
                  AA = (A( L )).toDouble();
                  // L(k-1,k-1) at A(k-1,k-1)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
               }
            }
         } else {
            // n is even
            if ( IFM == 1 ) {
               // A is normal
               if ( ILU == 0 ) {
                  // A is upper
                  for (J = 0; J <= K - 2; J++) {
                     zlassq(K-J-1, A( K+J+2+J*LDA ), 1, SCALE, S );
                  // L at A(k+1,0)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     zlassq(K+J, A( 0+J*LDA ), 1, SCALE, S );
                  // trap U at A(0,0)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = K;
                  // -> U(k,k) at A(k,0)
                  for (I = 0; I <= K - 1; I++) {
                     AA = (A( L )).toDouble();
                     // U(k+i,k+i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // U(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
               } else {
                  // ilu=1 & A is lower
                  for (J = 0; J <= K - 1; J++) {
                     zlassq(N-J-1, A( J+2+J*LDA ), 1, SCALE, S );
                     // trap L at A(1,0)
                  }
                  for (J = 1; J <= K - 1; J++) {
                     zlassq(J, A( 0+J*LDA ), 1, SCALE, S );
                     // U at A(0,0)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = 0;
                  // -> L(k,k) at A(0,0)
                  for (I = 0; I <= K - 1; I++) {
                     AA = (A( L )).toDouble();
                     // L(k-1+i,k-1+i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // L(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
               }
            } else {
               // A is xpose
               if ( ILU == 0 ) {
                  // A**H is upper
                  for (J = 1; J <= K - 1; J++) {
                     zlassq(J, A( 0+( K+1+J )*LDA ), 1, SCALE, S );
                  // U at A(0,k+1)
                  }
                  for (J = 0; J <= K - 1; J++) {
                     zlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                  // k by k rect. at A(0,0)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     zlassq(K-J-1, A( J+1+( J+K )*LDA ), 1, SCALE, S );
                  // L at A(0,k)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = 0 + K*LDA;
                  // -> U(k,k) at A(0,k)
                  AA = (A( L )).toDouble();
                  // U(k,k)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
                  L = L + LDA;
                  // -> U(0,0) at A(0,k+1)
                  for (J = K + 1; J <= N - 1; J++) {
                     AA = (A( L )).toDouble();
                     // -> U(j-k-1,j-k-1)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // -> U(j,j)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
                  // L=k-1+n*lda
                  // -> U(k-1,k-1) at A(k-1,n)
                  AA = (A( L )).toDouble();
                  // U(k,k)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
               } else {
                  // A**H is lower
                  for (J = 1; J <= K - 1; J++) {
                     zlassq(J, A( 0+( J+1 )*LDA ), 1, SCALE, S );
                  // U at A(0,1)
                  }
                  for (J = K + 1; J <= N; J++) {
                     zlassq(K, A( 0+J*LDA ), 1, SCALE, S );
                  // k by k rect. at A(0,k+1)
                  }
                  for (J = 0; J <= K - 2; J++) {
                     zlassq(K-J-1, A( J+1+J*LDA ), 1, SCALE, S );
                  // L at A(0,0)
                  }
                  S = S + S;
                  // double s for the off diagonal elements
                  L = 0;
                  // -> L(k,k) at A(0,0)
                  AA = (A( L )).toDouble();
                  // L(k,k) at A(0,0)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
                  L = LDA;
                  // -> L(0,0) at A(0,1)
                  for (I = 0; I <= K - 2; I++) {
                     AA = (A( L )).toDouble();
                     // L(i,i)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     AA = (A( L+1 )).toDouble();
                     // L(k+i+1,k+i+1)
                     if ( AA != ZERO ) {
                        if ( SCALE < AA ) {
                           S = ONE + S*( SCALE / AA )**2;
                           SCALE = AA;
                        } else {
                           S = S + ( AA / SCALE )**2;
                        }
                     }
                     L = L + LDA + 1;
                  }
                  // L-> k - 1 + k*lda or L(k-1,k-1) at A(k-1,k)
                  AA = (A( L )).toDouble();
                  // L(k-1,k-1) at A(k-1,k)
                  if ( AA != ZERO ) {
                     if ( SCALE < AA ) {
                        S = ONE + S*( SCALE / AA )**2;
                        SCALE = AA;
                     } else {
                        S = S + ( AA / SCALE )**2;
                     }
                  }
               }
            }
         }
         VALUE = SCALE*sqrt( S );
      }

      ZLANHF = VALUE;
      }
