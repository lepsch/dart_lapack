      SUBROUTINE CTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, LDB, LDT, N, M, L;
      // ..
      // .. Array Arguments ..
      COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX  ONE, ZERO
      const    ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int       I, J, P, MP, NP;
      COMPLEX   ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFG, CGEMV, CGERC, CTRMV, XERBLA
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
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -7
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('CTPQRT2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN

      for (I = 1; I <= N; I++) {

         // Generate elementary reflector H(I) to annihilate B(:,I)

         P = M-L+MIN( L, I )
         clarfg(P+1, A( I, I ), B( 1, I ), 1, T( I, 1 ) );
         if ( I.LT.N ) {

            // W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)]

            DO J = 1, N-I
               T( J, N ) = CONJG(A( I, I+J ))
            END DO
            cgemv('C', P, N-I, ONE, B( 1, I+1 ), LDB, B( 1, I ), 1, ONE, T( 1, N ), 1 );

            // C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H

            ALPHA = -CONJG(T( I, 1 ))
            DO J = 1, N-I
               A( I, I+J ) = A( I, I+J ) + ALPHA*CONJG(T( J, N ))
            END DO
            cgerc(P, N-I, ALPHA, B( 1, I ), 1, T( 1, N ), 1, B( 1, I+1 ), LDB );
         }
      END DO

      for (I = 2; I <= N; I++) {

         // T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I))

         ALPHA = -T( I, 1 )

         DO J = 1, I-1
            T( J, I ) = ZERO
         END DO
         P = MIN( I-1, L )
         MP = MIN( M-L+1, M )
         NP = MIN( P+1, N )

         // Triangular part of B2

         for (J = 1; J <= P; J++) {
            T( J, I ) = ALPHA*B( M-L+J, I )
         END DO
         ctrmv('U', 'C', 'N', P, B( MP, 1 ), LDB, T( 1, I ), 1 );

         // Rectangular part of B2

         cgemv('C', L, I-1-P, ALPHA, B( MP, NP ), LDB, B( MP, I ), 1, ZERO, T( NP, I ), 1 );

         // B1

         cgemv('C', M-L, I-1, ALPHA, B, LDB, B( 1, I ), 1, ONE, T( 1, I ), 1 );

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

         ctrmv('U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 );

         // T(I,I) = tau(I)

         T( I, I ) = T( I, 1 )
         T( I, 1 ) = ZERO
      END DO


      // End of CTPQRT2

      }
