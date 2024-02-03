      void stpttf(TRANSR, UPLO, N, AP, ARF, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( 0: * ), ARF( 0: * );

// =====================================================================

      // .. Parameters ..
      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K, NT;
      int                I, J, IJ;
      int                IJP, JP, LDA, JS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
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
      }
      if ( INFO != 0 ) {
         xerbla('STPTTF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( NORMALTRANSR ) {
            ARF( 0 ) = AP( 0 );
         } else {
            ARF( 0 ) = AP( 0 );
         }
         return;
      }

      // Size of array ARF(0:NT-1)

      NT = N*( N+1 ) / 2;

      // Set N1 and N2 depending on LOWER

      if ( LOWER ) {
         N2 = N / 2;
         N1 = N - N2;
      } else {
         N1 = N / 2;
         N2 = N - N1;
      }

      // If N is odd, set NISODD = true;
      // If N is even, set K = N/2 and NISODD = false;

      // set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe)
      // where noe = 0 if n is even, noe = 1 if n is odd

      if ( MOD( N, 2 ) == 0 ) {
         K = N / 2;
         NISODD = false;
         LDA = N + 1;
      } else {
         NISODD = true;
         LDA = N;
      }

      // ARF^C has lda rows and n+1-noe cols

      if ( !NORMALTRANSR) LDA = ( N+1 ) / 2;

      // start execution: there are eight cases

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               IJP = 0;
               JP = 0;
               for (J = 0; J <= N2; J++) {
                  for (I = J; I <= N - 1; I++) {
                     IJ = I + JP;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JP = JP + LDA;
               }
               for (I = 0; I <= N2 - 1; I++) {
                  for (J = 1 + I; J <= N2; J++) {
                     IJ = I + J*LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               IJP = 0;
               for (J = 0; J <= N1 - 1; J++) {
                  IJ = N2 + J;
                  for (I = 0; I <= J; I++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                     IJ = IJ + LDA;
                  }
               }
               JS = 0;
               for (J = N1; J <= N - 1; J++) {
                  IJ = JS;
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               IJP = 0;
               for (I = 0; I <= N2; I++) {
                  DO IJ = I*( LDA+1 ), N*LDA - 1, LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
               }
               JS = 1;
               for (J = 0; J <= N2 - 1; J++) {
                  for (IJ = JS; IJ <= JS + N2 - J - 1; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA + 1;
               }

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               IJP = 0;
               JS = N2*LDA;
               for (J = 0; J <= N1 - 1; J++) {
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }
               for (I = 0; I <= N1; I++) {
                  DO IJ = I, I + ( N1+I )*LDA, LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
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

               IJP = 0;
               JP = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (I = J; I <= N - 1; I++) {
                     IJ = 1 + I + JP;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JP = JP + LDA;
               }
               for (I = 0; I <= K - 1; I++) {
                  for (J = I; J <= K - 1; J++) {
                     IJ = I + J*LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               IJP = 0;
               for (J = 0; J <= K - 1; J++) {
                  IJ = K + 1 + J;
                  for (I = 0; I <= J; I++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                     IJ = IJ + LDA;
                  }
               }
               JS = 0;
               for (J = K; J <= N - 1; J++) {
                  IJ = JS;
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }

            }

         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               IJP = 0;
               for (I = 0; I <= K - 1; I++) {
                  DO IJ = I + ( I+1 )*LDA, ( N+1 )*LDA - 1, LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
               }
               JS = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (IJ = JS; IJ <= JS + K - J - 1; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA + 1;
               }

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               IJP = 0;
               JS = ( K+1 )*LDA;
               for (J = 0; J <= K - 1; J++) {
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }
               for (I = 0; I <= K - 1; I++) {
                  DO IJ = I, I + ( K+I )*LDA, LDA;
                     ARF( IJ ) = AP( IJP );
                     IJP = IJP + 1;
                  }
               }

            }

         }

      }

      return;
      }
