      void ctrttf(TRANSR, UPLO, N, A, LDA, ARF, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      Complex            A( 0: LDA-1, 0: * ), ARF( 0: * );
      // ..

// =====================================================================

      // .. Parameters ..
      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                I, IJ, J, K, L, N1, N2, NT, NX2, NP1X2;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CTRTTF', -INFO );
         return;
      }

      // Quick return if possible

      if ( N <= 1 ) {
         if ( N == 1 ) {
            if ( NORMALTRANSR ) {
               ARF[0] = A( 0, 0 );
            } else {
               ARF[0] = CONJG( A( 0, 0 ) );
            }
         }
         return;
      }

      // Size of array ARF(1:2,0:nt-1)

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

              // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
              // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
              // T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n

               IJ = 0;
               for (J = 0; J <= N2; J++) {
                  for (I = N1; I <= N2 + J; I++) {
                     ARF[IJ] = CONJG( A( N2+J, I ) );
                     IJ = IJ + 1;
                  }
                  for (I = J; I <= N - 1; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
               }

            } else {

              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n

               IJ = NT - N;
               for (J = N - 1; J >= N1; J--) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
                  for (L = J - N1; L <= N1 - 1; L++) {
                     ARF[IJ] = CONJG( A( J-N1, L ) );
                     IJ = IJ + 1;
                  }
                  IJ = IJ - NX2;
               }

            }

         } else {

            // N is odd and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
               // T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1

               IJ = 0;
               for (J = 0; J <= N2 - 1; J++) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
                  for (I = N1 + J; I <= N - 1; I++) {
                     ARF[IJ] = A( I, N1+J );
                     IJ = IJ + 1;
                  }
               }
               for (J = N2; J <= N - 1; J++) {
                  for (I = 0; I <= N1 - 1; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
               // T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda=n2

               IJ = 0;
               for (J = 0; J <= N1; J++) {
                  for (I = N1; I <= N - 1; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
               }
               for (J = 0; J <= N1 - 1; J++) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
                  for (L = N2 + J; L <= N - 1; L++) {
                     ARF[IJ] = CONJG( A( N2+J, L ) );
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

               // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
               // T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1

               IJ = 0;
               for (J = 0; J <= K - 1; J++) {
                  for (I = K; I <= K + J; I++) {
                     ARF[IJ] = CONJG( A( K+J, I ) );
                     IJ = IJ + 1;
                  }
                  for (I = J; I <= N - 1; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1

               IJ = NT - N - 1;
               for (J = N - 1; J >= K; J--) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
                  for (L = J - K; L <= K - 1; L++) {
                     ARF[IJ] = CONJG( A( J-K, L ) );
                     IJ = IJ + 1;
                  }
                  IJ = IJ - NP1X2;
               }

            }

         } else {

            // N is even and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B)
               // T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) :
               // T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k

               IJ = 0;
               J = K;
               for (I = K; I <= N - 1; I++) {
                  ARF[IJ] = A( I, J );
                  IJ = IJ + 1;
               }
               for (J = 0; J <= K - 2; J++) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
                  for (I = K + 1 + J; I <= N - 1; I++) {
                     ARF[IJ] = A( I, K+1+J );
                     IJ = IJ + 1;
                  }
               }
               for (J = K - 1; J <= N - 1; J++) {
                  for (I = 0; I <= K - 1; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
               }

            } else {

               // SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B)
               // T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0)
               // T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k

               IJ = 0;
               for (J = 0; J <= K; J++) {
                  for (I = K; I <= N - 1; I++) {
                     ARF[IJ] = CONJG( A( J, I ) );
                     IJ = IJ + 1;
                  }
               }
               for (J = 0; J <= K - 2; J++) {
                  for (I = 0; I <= J; I++) {
                     ARF[IJ] = A( I, J );
                     IJ = IJ + 1;
                  }
                  for (L = K + 1 + J; L <= N - 1; L++) {
                     ARF[IJ] = CONJG( A( K+1+J, L ) );
                     IJ = IJ + 1;
                  }
               }

               // Note that here J = K-1

               for (I = 0; I <= J; I++) {
                  ARF[IJ] = A( I, J );
                  IJ = IJ + 1;
               }

            }

         }

      }

      return;
      }