      SUBROUTINE CGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDT, M, N, MB;
      // ..
      // .. Array Arguments ..
      COMPLEX A( LDA, * ), T( LDT, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, IINFO, K;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQT3, CLARFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( MB.LT.1 .OR. (MB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 )) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDT.LT.MB ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         xerbla('CGELQT', -INFO );
         RETURN
      }

      // Quick return if possible

      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN

      // Blocked loop of length K

      DO I = 1, K,  MB
         IB = MIN( K-I+1, MB )

      // Compute the LQ factorization of the current block A(I:M,I:I+IB-1)

         cgelqt3(IB, N-I+1, A(I,I), LDA, T(1,I), LDT, IINFO );
         if ( I+IB.LE.M ) {

      // Update by applying H**T to A(I:M,I+IB:N) from the right

         clarfb('R', 'N', 'F', 'R', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, T( 1, I ), LDT, A( I+IB, I ), LDA, WORK , M-I-IB+1 );
         }
      END DO
      RETURN

      // End of CGELQT

      }
