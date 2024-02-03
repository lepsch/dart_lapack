      SUBROUTINE CGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQ2, CLARFB, CLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      K = MIN( M, N )
      NB = ILAENV( 1, 'CGELQF', ' ', M, N, -1, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      } else if ( .NOT.LQUERY ) {
         IF( LWORK.LE.0 .OR. ( N.GT.0 .AND. LWORK.LT.MAX( 1, M ) ) ) INFO = -7
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGELQF', -INFO )
         RETURN
      } else if ( LQUERY ) {
         if ( K.EQ.0 ) {
            LWKOPT = 1
         } else {
            LWKOPT = M*NB
         }
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      }

      // Quick return if possible

      if ( K.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NBMIN = 2
      NX = 0
      IWS = M
      if ( NB.GT.1 .AND. NB.LT.K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'CGELQF', ' ', M, N, -1, -1 ) )
         if ( NX.LT.K ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M
            IWS = LDWORK*NB
            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'CGELQF', ' ', M, N, -1, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) {

         // Use blocked code initially

         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )

            // Compute the LQ factorization of the current block
            // A(i:i+ib-1,i:n)

            CALL CGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )
            if ( I+IB.LE.M ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               CALL CLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, LDWORK )

               // Apply H to A(i+ib:m,i:n) from the right

               CALL CLARFB( 'Right', 'No transpose', 'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ), LDWORK )
            }
   10    CONTINUE
      } else {
         I = 1
      }

      // Use unblocked code to factor the last or only block.

      IF( I.LE.K ) CALL CGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

      WORK( 1 ) = SROUNDUP_LWORK( IWS )
      RETURN

      // End of CGELQF

      }
