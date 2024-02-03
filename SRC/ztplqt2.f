      SUBROUTINE ZTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int            INFO, LDA, LDB, LDT, N, M, L;
      // ..
      // .. Array Arguments ..
      COMPLEX*16     A( LDA, * ), B( LDB, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16  ONE, ZERO
      const    ZERO = ( 0.0, 0.0 ),ONE  = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int       I, J, P, MP, NP;
      COMPLEX*16   ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARFG, ZGEMV, ZGERC, ZTRMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( L < 0 || L > MIN(M,N) ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, M ) ) {
         INFO = -7
      } else if ( LDT < MAX( 1, M ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('ZTPLQT2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || M == 0) RETURN;

      for (I = 1; I <= M; I++) {

         // Generate elementary reflector H(I) to annihilate B(I,:)

         P = N-L+MIN( L, I )
         zlarfg(P+1, A( I, I ), B( I, 1 ), LDB, T( 1, I ) );
         T(1,I)=CONJG(T(1,I))
         if ( I < M ) {
            for (J = 1; J <= P; J++) {
               B( I, J ) = CONJG(B(I,J))
            }

            // W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]

            for (J = 1; J <= M-I; J++) {
               T( M, J ) = (A( I+J, I ))
            }
            zgemv('N', M-I, P, ONE, B( I+1, 1 ), LDB, B( I, 1 ), LDB, ONE, T( M, 1 ), LDT );

            // C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H

            ALPHA = -(T( 1, I ))
            for (J = 1; J <= M-I; J++) {
               A( I+J, I ) = A( I+J, I ) + ALPHA*(T( M, J ))
            }
            zgerc(M-I, P, (ALPHA),  T( M, 1 ), LDT, B( I, 1 ), LDB, B( I+1, 1 ), LDB );
            for (J = 1; J <= P; J++) {
               B( I, J ) = CONJG(B(I,J))
            }
         }
      }

      for (I = 2; I <= M; I++) {

         // T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N))

         ALPHA = -(T( 1, I ))
         for (J = 1; J <= I-1; J++) {
            T( I, J ) = ZERO
         }
         P = MIN( I-1, L )
         NP = MIN( N-L+1, N )
         MP = MIN( P+1, M )
         for (J = 1; J <= N-L+P; J++) {
           B(I,J)=CONJG(B(I,J))
         }

         // Triangular part of B2

         for (J = 1; J <= P; J++) {
            T( I, J ) = (ALPHA*B( I, N-L+J ))
         }
         ztrmv('L', 'N', 'N', P, B( 1, NP ), LDB, T( I, 1 ), LDT );

         // Rectangular part of B2

         zgemv('N', I-1-P, L,  ALPHA, B( MP, NP ), LDB, B( I, NP ), LDB, ZERO, T( I,MP ), LDT );

         // B1


         zgemv('N', I-1, N-L, ALPHA, B, LDB, B( I, 1 ), LDB, ONE, T( I, 1 ), LDT );



         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)

         for (J = 1; J <= I-1; J++) {
            T(I,J)=CONJG(T(I,J))
         }
         ztrmv('L', 'C', 'N', I-1, T, LDT, T( I, 1 ), LDT );
         for (J = 1; J <= I-1; J++) {
            T(I,J)=CONJG(T(I,J))
         }
         for (J = 1; J <= N-L+P; J++) {
            B(I,J)=CONJG(B(I,J))
         }

         // T(I,I) = tau(I)

         T( I, I ) = T( 1, I )
         T( 1, I ) = ZERO
      }
      for (I = 1; I <= M; I++) {
         for (J = I+1; J <= M; J++) {
            T(I,J)=(T(J,I))
            T(J,I)=ZERO
         }
      }


      // End of ZTPLQT2

      }
