      void ctpttf(TRANSR, UPLO, N, AP, ARF, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( 0: * ), ARF( 0: * );

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
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = LSAME( TRANSR, 'N' );
      LOWER = LSAME( UPLO, 'L' );
      if ( !NORMALTRANSR && !LSAME( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('CTPTTF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( NORMALTRANSR ) {
            ARF[0] = AP( 0 );
         } else {
            ARF[0] = CONJG( AP( 0 ) );
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

      if ( (N % 2) == 0 ) {
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

              // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
              // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
              // T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n

               IJP = 0;
               JP = 0;
               for (J = 0; J <= N2; J++) {
                  for (I = J; I <= N - 1; I++) {
                     IJ = I + JP;
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JP = JP + LDA;
               }
               for (I = 0; I <= N2 - 1; I++) {
                  for (J = 1 + I; J <= N2; J++) {
                     IJ = I + J*LDA;
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                  }
               }

            } else {

              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)

               IJP = 0;
               for (J = 0; J <= N1 - 1; J++) {
                  IJ = N2 + J;
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                     IJ = IJ + LDA;
                  }
               }
               JS = 0;
               for (J = N1; J <= N - 1; J++) {
                  IJ = JS;
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }

            }

         } else {

            // N is odd and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
               // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1

               IJP = 0;
               for (I = 0; I <= N2; I++) {
                  for (IJ = I*( LDA+1 ); LDA < 0 ? IJ >= N*LDA - 1 : IJ <= N*LDA - 1; IJ += LDA) {
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                  }
               }
               JS = 1;
               for (J = 0; J <= N2 - 1; J++) {
                  for (IJ = JS; IJ <= JS + N2 - J - 1; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA + 1;
               }

            } else {

               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
               // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2

               IJP = 0;
               JS = N2*LDA;
               for (J = 0; J <= N1 - 1; J++) {
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }
               for (I = 0; I <= N1; I++) {
                  for (IJ = I; LDA < 0 ? IJ >= I + ( N1+I )*LDA : IJ <= I + ( N1+I )*LDA; IJ += LDA) {
                     ARF[IJ] = CONJG( AP( IJP ) );
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

               // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
               // T1 -> a(1), T2 -> a(0), S -> a(k+1)

               IJP = 0;
               JP = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (I = J; I <= N - 1; I++) {
                     IJ = 1 + I + JP;
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JP = JP + LDA;
               }
               for (I = 0; I <= K - 1; I++) {
                  for (J = I; J <= K - 1; J++) {
                     IJ = I + J*LDA;
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                  }
               }

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)

               IJP = 0;
               for (J = 0; J <= K - 1; J++) {
                  IJ = K + 1 + J;
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                     IJ = IJ + LDA;
                  }
               }
               JS = 0;
               for (J = K; J <= N - 1; J++) {
                  IJ = JS;
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }

            }

         } else {

            // N is even and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

               IJP = 0;
               for (I = 0; I <= K - 1; I++) {
                  DO IJ = I + ( I+1 )*LDA, ( N+1 )*LDA - 1, LDA;
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                  }
               }
               JS = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (IJ = JS; IJ <= JS + K - J - 1; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA + 1;
               }

            } else {

               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

               IJP = 0;
               JS = ( K+1 )*LDA;
               for (J = 0; J <= K - 1; J++) {
                  for (IJ = JS; IJ <= JS + J; IJ++) {
                     ARF[IJ] = AP( IJP );
                     IJP = IJP + 1;
                  }
                  JS = JS + LDA;
               }
               for (I = 0; I <= K - 1; I++) {
                  for (IJ = I; LDA < 0 ? IJ >= I + ( K+I )*LDA : IJ <= I + ( K+I )*LDA; IJ += LDA) {
                     ARF[IJ] = CONJG( AP( IJP ) );
                     IJP = IJP + 1;
                  }
               }

            }

         }

      }

      return;
      }
