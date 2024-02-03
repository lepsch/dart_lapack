      void stfttr(TRANSR, UPLO, N, ARF, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      REAL               A( 0: LDA-1, 0: * ), ARF( 0: * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K, NT, NX2, NP1X2;
      int                I, J, L, IJ;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = LSAME( TRANSR, 'N' );
      LOWER = LSAME( UPLO, 'L' );
      if ( !NORMALTRANSR && !LSAME( TRANSR, 'T' ) ) {
         INFO = -1;
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('STFTTR', -INFO );
         return;
      }

      // Quick return if possible

      if ( N <= 1 ) {
         if ( N == 1 ) {
            A( 0, 0 ) = ARF( 0 );
         }
         return;
      }

      // Size of array ARF(0:nt-1)

      NT = N*( N+1 ) / 2;

      // set N1 and N2 depending on LOWER: for N even N1=N2=K

      if ( LOWER ) {
         N2 = N / 2;
         N1 = N - N2;
      } else {
         N1 = N / 2;
         N2 = N - N1;
      }

      // If N is odd, set NISODD = true , LDA=N+1 and A is (N+1)--by--K2.
      // If N is even, set K = N/2 and NISODD = false , LDA=N and A is
      // N--by--(N+1)/2.

      if ( (N % 2) == 0 ) {
         K = N / 2;
         NISODD = false;
         if ( !LOWER) NP1X2 = N + N + 2;
      } else {
         NISODD = true;
         if ( !LOWER) NX2 = N + N;
      }

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               IJ = 0;
               for (J = 0; J <= N2; J++) {
                  for (I = N1; I <= N2 + J; I++) {
                     A( N2+J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (I = J; I <= N - 1; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N;
               for (J = N - 1; J >= N1; J--) { //
                  for (I = 0; I <= J; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (L = J - N1; L <= N1 - 1; L++) {
                     A( J-N1, L ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  IJ = IJ - NX2;
               }

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               IJ = 0;
               for (J = 0; J <= N2 - 1; J++) {
                  for (I = 0; I <= J; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (I = N1 + J; I <= N - 1; I++) {
                     A( I, N1+J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }
               for (J = N2; J <= N - 1; J++) {
                  for (I = 0; I <= N1 - 1; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               IJ = 0;
               for (J = 0; J <= N1; J++) {
                  for (I = N1; I <= N - 1; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }
               for (J = 0; J <= N1 - 1; J++) {
                  for (I = 0; I <= J; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (L = N2 + J; L <= N - 1; L++) {
                     A( N2+J, L ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }

            }

         }

      } else {

         // N is even

         if ( NORMALTRANSR ) {

            // N is even and TRANSR = 'N'

            if ( LOWER ) {

               // N is even, TRANSR = 'N', and UPLO = 'L'

               IJ = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (I = K; I <= K + J; I++) {
                     A( K+J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (I = J; I <= N - 1; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N - 1;
               for (J = N - 1; J >= K; J--) { //
                  for (I = 0; I <= J; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (L = J - K; L <= K - 1; L++) {
                     A( J-K, L ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  IJ = IJ - NP1X2;
               }

            }

         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               IJ = 0;
               J = K;
               for (I = K; I <= N - 1; I++) {
                  A( I, J ) = ARF( IJ );
                  IJ = IJ + 1;
               }
               for (J = 0; J <= K - 2; J++) {
                  for (I = 0; I <= J; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (I = K + 1 + J; I <= N - 1; I++) {
                     A( I, K+1+J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }
               for (J = K - 1; J <= N - 1; J++) {
                  for (I = 0; I <= K - 1; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               IJ = 0;
               for (J = 0; J <= K; J++) {
                  for (I = K; I <= N - 1; I++) {
                     A( J, I ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }
               for (J = 0; J <= K - 2; J++) {
                  for (I = 0; I <= J; I++) {
                     A( I, J ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
                  for (L = K + 1 + J; L <= N - 1; L++) {
                     A( K+1+J, L ) = ARF( IJ );
                     IJ = IJ + 1;
                  }
               }
               // Note that here, on exit of the loop, J = K-1
               for (I = 0; I <= J; I++) {
                  A( I, J ) = ARF( IJ );
                  IJ = IJ + 1;
               }

            }

         }

      }

      return;
      }
