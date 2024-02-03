      SUBROUTINE CGEQRF ( M, N, A, LDA, TAU, WORK, LWORK, INFO )

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
      int                I, IB, IINFO, IWS, J, K, LWKOPT, NB, NBMIN, NX, LBWORK, NT, LLWORK;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQR2, CLARFB, CLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, MAX, MIN, REAL
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      INFO = 0
      NBMIN = 2
      NX = 0
      IWS = N
      K = MIN( M, N )
      NB = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )

      if ( NB.GT.1 .AND. NB.LT.K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'CGEQRF', ' ', M, N, -1, -1 ) )
      }

      // Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.:

             // NB=3     2NB=6       K=10
             // |        |           |
       // 1--2--3--4--5--6--7--8--9--10
                   // |     \________/
                // K-NX=5      NT=4

      // So here 4 x 4 is the last T stored in the workspace

      NT = K-CEILING(REAL(K-NX)/REAL(NB))*NB


      // optimal workspace = space for dlarfb + space for normal T's + space for the last T

      LLWORK = MAX (MAX((N-M)*K, (N-M)*NB), MAX(K*NB, NB*NB))
      LLWORK = CEILING(REAL(LLWORK)/REAL(NB))

      if ( K.EQ.0 ) {

         LBWORK = 0
         LWKOPT = 1
         WORK( 1 ) = LWKOPT

      } else if ( NT.GT.NB ) {

          LBWORK = K-NT

          // Optimal workspace for dlarfb = MAX(1,N)*NT

          LWKOPT = (LBWORK+LLWORK)*NB
          WORK( 1 ) = SROUNDUP_LWORK(LWKOPT+NT*NT)

      } else {

          LBWORK = CEILING(REAL(K)/REAL(NB))*NB
          LWKOPT = (LBWORK+LLWORK-NB)*NB
          WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

      }


      // Test the input arguments

      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      } else if ( .NOT.LQUERY ) {
         IF( LWORK.LE.0 .OR. ( M.GT.0 .AND. LWORK.LT.MAX( 1, N ) ) ) INFO = -7
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGEQRF', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( K.EQ.0 ) {
         RETURN
      }

      if ( NB.GT.1 .AND. NB.LT.K ) {

         if ( NX.LT.K ) {

            // Determine if workspace is large enough for blocked code.

            if ( NT.LE.NB ) {
                IWS = (LBWORK+LLWORK-NB)*NB
            } else {
                IWS = (LBWORK+LLWORK)*NB+NT*NT
            }

            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               if ( NT.LE.NB ) {
                    NB = LWORK / (LLWORK+(LBWORK-NB))
               } else {
                    NB = (LWORK-NT*NT)/(LBWORK+LLWORK)
               }
                NBMIN = MAX( 2, ILAENV( 2, 'CGEQRF', ' ', M, N, -1, -1 ) )
            }
         }
      }

      if ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) {

         // Use blocked code initially

         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )

            // Update the current column using old T's

            DO 20 J = 1, I - NB, NB

               // Apply H' to A(J:M,I:I+IB-1) from the left

               CALL CLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, IB, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), IB)

20          CONTINUE

            // Compute the QR factorization of the current block
            // A(I:M,I:I+IB-1)

            CALL CGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1), IINFO )

            if ( I+IB.LE.N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               CALL CLARFT( 'Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK )

            }
   10    CONTINUE
      } else {
         I = 1
      }

      // Use unblocked code to factor the last or only block.

      if ( I.LE.K ) {

         if ( I .NE. 1 ) {

             DO 30 J = 1, I - NB, NB

                 // Apply H' to A(J:M,I:K) from the left

                 CALL CLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, K-I+1, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), K-I+1)
30           CONTINUE
              CALL CGEQR2( M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1),IINFO )

         } else {

         // Use unblocked code to factor the last or only block.

         CALL CGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO )

         }
      }



      // Apply update to the column M+1:N when N > M

      if ( M.LT.N .AND. I.NE.1) {

          // Form the last triangular factor of the block reflector
          // H = H(i) H(i+1) . . . H(i+ib-1)

          if ( NT .LE. NB ) {
               CALL CLARFT( 'Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK )
          } else {
               CALL CLARFT( 'Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+1), NT )
          }


          // Apply H' to A(1:M,M+1:N) from the left

          DO 40 J = 1, K-NX, NB

               IB = MIN( K-J+1, NB )
                CALL CLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, IB, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)

40       CONTINUE

         if ( NT.LE.NB ) {
             CALL CLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)
         } else {
             CALL CLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(LBWORK*NB+1), NT, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M)
         }

      }

      WORK( 1 ) = SROUNDUP_LWORK(IWS)
      RETURN

      // End of CGEQRF

      }
