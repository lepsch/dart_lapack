      SUBROUTINE STPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, LDB, LDT, N, M, L;
      // ..
      // .. Array Arguments ..
      REAL   A( LDA, * ), B( LDB, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL  ONE, ZERO
      const    ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int       I, J, P, MP, NP;
      REAL   ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, SGEMV, SGER, STRMV, XERBLA
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
         xerbla('STPLQT2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN

      for (I = 1; I <= M; I++) {

         // Generate elementary reflector H(I) to annihilate B(I,:)

         P = N-L+MIN( L, I )
         slarfg(P+1, A( I, I ), B( I, 1 ), LDB, T( 1, I ) );
         if ( I.LT.M ) {

            // W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]

            DO J = 1, M-I
               T( M, J ) = (A( I+J, I ))
            END DO
            sgemv('N', M-I, P, ONE, B( I+1, 1 ), LDB, B( I, 1 ), LDB, ONE, T( M, 1 ), LDT );

            // C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H

            ALPHA = -(T( 1, I ))
            DO J = 1, M-I
               A( I+J, I ) = A( I+J, I ) + ALPHA*(T( M, J ))
            END DO
            sger(M-I, P, ALPHA,  T( M, 1 ), LDT, B( I, 1 ), LDB, B( I+1, 1 ), LDB );
         }
      END DO

      for (I = 2; I <= M; I++) {

         // T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H)

         ALPHA = -T( 1, I )

         DO J = 1, I-1
            T( I, J ) = ZERO
         END DO
         P = MIN( I-1, L )
         NP = MIN( N-L+1, N )
         MP = MIN( P+1, M )

         // Triangular part of B2

         for (J = 1; J <= P; J++) {
            T( I, J ) = ALPHA*B( I, N-L+J )
         END DO
         strmv('L', 'N', 'N', P, B( 1, NP ), LDB, T( I, 1 ), LDT );

         // Rectangular part of B2

         sgemv('N', I-1-P, L,  ALPHA, B( MP, NP ), LDB, B( I, NP ), LDB, ZERO, T( I,MP ), LDT );

         // B1

         sgemv('N', I-1, N-L, ALPHA, B, LDB, B( I, 1 ), LDB, ONE, T( I, 1 ), LDT );

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)

        strmv('L', 'T', 'N', I-1, T, LDT, T( I, 1 ), LDT );

         // T(I,I) = tau(I)

         T( I, I ) = T( 1, I )
         T( 1, I ) = ZERO
      END DO
      for (I = 1; I <= M; I++) {
         DO J= I+1,M
            T(I,J)=T(J,I)
            T(J,I)= ZERO
         END DO
      END DO


      // End of STPLQT2

      }
