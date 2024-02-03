      SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQR2, ZLARFB, ZLARFT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      K = MIN( M, N )
      INFO = 0
      NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4
      } else if ( .NOT.LQUERY ) {
         IF( LWORK.LE.0 || ( M > 0 && LWORK < MAX( 1, N ) ) ) INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('ZGEQRF', -INFO );
         RETURN
      } else if ( LQUERY ) {
         if ( K == 0 ) {
            LWKOPT = 1
         } else {
            LWKOPT = N*NB
         }
         WORK( 1 ) = LWKOPT
         RETURN
      }

      // Quick return if possible

      if ( K == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NBMIN = 2
      NX = 0
      IWS = N
      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEQRF', ' ', M, N, -1, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN && NB < K && NX < K ) {

         // Use blocked code initially

         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )

            // Compute the QR factorization of the current block
            // A(i:m,i:i+ib-1)

            zgeqr2(M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO );
            if ( I+IB.LE.N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               zlarft('Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H**H to A(i:m,i+ib:n) from the left

               zlarfb('Left', 'Conjugate transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK );
            }
         } // 10
      } else {
         I = 1
      }

      // Use unblocked code to factor the last or only block.

      if (I.LE.K) CALL ZGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO );

      WORK( 1 ) = IWS
      RETURN

      // End of ZGEQRF

      }
