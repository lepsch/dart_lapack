      SUBROUTINE SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDT, M, N, NB;
      // ..
      // .. Array Arguments ..
      REAL A( LDA, * ), T( LDT, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, IINFO, K;
      bool       USE_RECURSIVE_QR;
      const    USE_RECURSIVE_QR= true ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRT2, SGEQRT3, SLARFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NB.LT.1 || ( NB.GT.MIN(M,N) && MIN(M,N).GT.0 ) ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDT.LT.NB ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('SGEQRT', -INFO );
         RETURN
      }

      // Quick return if possible

      K = MIN( M, N )
      if (K == 0) RETURN;

      // Blocked loop of length K

      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )

      // Compute the QR factorization of the current block A(I:M,I:I+IB-1)

         if ( USE_RECURSIVE_QR ) {
            sgeqrt3(M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO );
         } else {
            sgeqrt2(M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO );
         }
         if ( I+IB.LE.N ) {

      // Update by applying H**T to A(I:M,I+IB:N) from the left

            slarfb('L', 'T', 'F', 'C', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, T( 1, I ), LDT, A( I, I+IB ), LDA, WORK , N-I-IB+1 );
         }
      }
      RETURN

      // End of SGEQRT

      }
