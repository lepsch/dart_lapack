      void ztpqrt2(M, N, L, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> T, final int LDT, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int       INFO, LDA, LDB, LDT, N, M, L;
      Complex   A( LDA, * ), B( LDB, * ), T( LDT, * );
      // ..

      Complex  ONE, ZERO;
      const    ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      int       I, J, P, MP, NP;
      Complex   ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARFG, ZGEMV, ZGERC, ZTRMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( L < 0 || L > min(M,N) ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -7;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('ZTPQRT2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || M == 0) return;

      for (I = 1; I <= N; I++) {

         // Generate elementary reflector H(I) to annihilate B(:,I)

         P = M-L+min( L, I );
         zlarfg(P+1, A( I, I ), B( 1, I ), 1, T( I, 1 ) );
         if ( I < N ) {

            // W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)]

            for (J = 1; J <= N-I; J++) {
               T[J][N] = CONJG(A( I, I+J ));
            }
            zgemv('C', P, N-I, ONE, B( 1, I+1 ), LDB, B( 1, I ), 1, ONE, T( 1, N ), 1 );

            // C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H

            ALPHA = -CONJG(T( I, 1 ));
            for (J = 1; J <= N-I; J++) {
               A[I][I+J] = A( I, I+J ) + ALPHA*CONJG(T( J, N ));
            }
            zgerc(P, N-I, ALPHA, B( 1, I ), 1, T( 1, N ), 1, B( 1, I+1 ), LDB );
         }
      }

      for (I = 2; I <= N; I++) {

         // T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I))

         ALPHA = -T( I, 1 );

         for (J = 1; J <= I-1; J++) {
            T[J][I] = ZERO;
         }
         P = min( I-1, L );
         MP = min( M-L+1, M );
         NP = min( P+1, N );

         // Triangular part of B2

         for (J = 1; J <= P; J++) {
            T[J][I] = ALPHA*B( M-L+J, I );
         }
         ztrmv('U', 'C', 'N', P, B( MP, 1 ), LDB, T( 1, I ), 1 );

         // Rectangular part of B2

         zgemv('C', L, I-1-P, ALPHA, B( MP, NP ), LDB, B( MP, I ), 1, ZERO, T( NP, I ), 1 );

         // B1

         zgemv('C', M-L, I-1, ALPHA, B, LDB, B( 1, I ), 1, ONE, T( 1, I ), 1 );

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

         ztrmv('U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 );

         // T(I,I) = tau(I)

         T[I][I] = T( I, 1 );
         T[I][1] = ZERO;
      }

      }
