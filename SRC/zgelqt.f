      SUBROUTINE ZGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDT, M, N, MB;
      // ..
      // .. Array Arguments ..
      COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, IINFO, K;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGELQT3, ZLARFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( MB < 1 || (MB > MIN(M,N) && MIN(M,N) > 0 )) {
         INFO = -3;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5;
      } else if ( LDT < MB ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('ZGELQT', -INFO );
         return;
      }

      // Quick return if possible

      K = MIN( M, N );
      if (K == 0) RETURN;

      // Blocked loop of length K

      DO I = 1, K,  MB;
         IB = MIN( K-I+1, MB );

      // Compute the LQ factorization of the current block A(I:M,I:I+IB-1)

         zgelqt3(IB, N-I+1, A(I,I), LDA, T(1,I), LDT, IINFO );
         if ( I+IB <= M ) {

      // Update by applying H**T to A(I:M,I+IB:N) from the right

         zlarfb('R', 'N', 'F', 'R', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, T( 1, I ), LDT, A( I+IB, I ), LDA, WORK , M-I-IB+1 );
         }
      }
      return;

      // End of ZGELQT

      }
