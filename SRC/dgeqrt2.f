      SUBROUTINE DGEQRT2( M, N, A, LDA, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, LDT, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), T( LDT, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double            ONE, ZERO;
      const    ONE = 1.0D+00, ZERO = 0.0D+00 ;
      // ..
      // .. Local Scalars ..
      int       I, K;
      double             AII, ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DGEMV, DGER, DTRMV, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -2
      } else if ( M.LT.N ) {
         INFO = -1
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DGEQRT2', -INFO )
         RETURN
      }

      K = MIN( M, N )

      DO I = 1, K

         // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)

         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, T( I, 1 ) )
         if ( I.LT.N ) {

            // Apply H(i) to A(I:M,I+1:N) from the left

            AII = A( I, I )
            A( I, I ) = ONE

            // W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]

            CALL DGEMV( 'T',M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, T( 1, N ), 1 )

            // A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H

            ALPHA = -(T( I, 1 ))
            CALL DGER( M-I+1, N-I, ALPHA, A( I, I ), 1, T( 1, N ), 1, A( I, I+1 ), LDA )
            A( I, I ) = AII
         }
      END DO

      DO I = 2, N
         AII = A( I, I )
         A( I, I ) = ONE

         // T(1:I-1,I) := alpha * A(I:M,1:I-1)**T * A(I:M,I)

         ALPHA = -T( I, 1 )
         CALL DGEMV( 'T', M-I+1, I-1, ALPHA, A( I, 1 ), LDA, A( I, I ), 1, ZERO, T( 1, I ), 1 )
         A( I, I ) = AII

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

         CALL DTRMV( 'U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 )

            // T(I,I) = tau(I)

            T( I, I ) = T( I, 1 )
            T( I, 1) = ZERO
      END DO


      // End of DGEQRT2

      }
