      SUBROUTINE DTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, LDB, LDT, N, M, L;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), T( LDT, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double            ONE, ZERO;
      const    ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int       I, J, P, MP, NP;
      double             ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DGEMV, DGER, DTRMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( L.LT.0 .OR. L.GT.MIN(M,N) ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -7
      } else if ( LDT.LT.MAX( 1, M ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('DTPLQT2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0 .OR. M.EQ.0) RETURN;

      for (I = 1; I <= M; I++) {

         // Generate elementary reflector H(I) to annihilate B(I,:)

         P = N-L+MIN( L, I )
         dlarfg(P+1, A( I, I ), B( I, 1 ), LDB, T( 1, I ) );
         if ( I.LT.M ) {

            // W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]

            for (J = 1; J <= M-I; J++) {
               T( M, J ) = (A( I+J, I ))
            }
            dgemv('N', M-I, P, ONE, B( I+1, 1 ), LDB, B( I, 1 ), LDB, ONE, T( M, 1 ), LDT );

            // C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H

            ALPHA = -(T( 1, I ))
            for (J = 1; J <= M-I; J++) {
               A( I+J, I ) = A( I+J, I ) + ALPHA*(T( M, J ))
            }
            dger(M-I, P, ALPHA,  T( M, 1 ), LDT, B( I, 1 ), LDB, B( I+1, 1 ), LDB );
         }
      }

      for (I = 2; I <= M; I++) {

         // T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H)

         ALPHA = -T( 1, I )

         for (J = 1; J <= I-1; J++) {
            T( I, J ) = ZERO
         }
         P = MIN( I-1, L )
         NP = MIN( N-L+1, N )
         MP = MIN( P+1, M )

         // Triangular part of B2

         for (J = 1; J <= P; J++) {
            T( I, J ) = ALPHA*B( I, N-L+J )
         }
         dtrmv('L', 'N', 'N', P, B( 1, NP ), LDB, T( I, 1 ), LDT );

         // Rectangular part of B2

         dgemv('N', I-1-P, L,  ALPHA, B( MP, NP ), LDB, B( I, NP ), LDB, ZERO, T( I,MP ), LDT );

         // B1

         dgemv('N', I-1, N-L, ALPHA, B, LDB, B( I, 1 ), LDB, ONE, T( I, 1 ), LDT );

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)

        dtrmv('L', 'T', 'N', I-1, T, LDT, T( I, 1 ), LDT );

         // T(I,I) = tau(I)

         T( I, I ) = T( 1, I )
         T( 1, I ) = ZERO
      }
      for (I = 1; I <= M; I++) {
         for (J = I+1; J <= M; J++) {
            T(I,J)=T(J,I)
            T(J,I)= ZERO
         }
      }


      // End of DTPLQT2

      }
