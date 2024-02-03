      SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IWS, KI, KK, LDWORK, LWKMIN, LWKOPT, M1, MU, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARZB, ZLARZT, ZLATRZ
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

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }

      if ( INFO == 0 ) {
         if ( M == 0 || M == N ) {
            LWKOPT = 1;
            LWKMIN = 1;
         } else {

            // Determine the block size.

            NB = ILAENV( 1, 'ZGERQF', ' ', M, N, -1, -1 );
            LWKOPT = M*NB;
            LWKMIN = MAX( 1, M );
         }
         WORK( 1 ) = LWKOPT;

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -7;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZTZRZF', -INFO );
         RETURN;
      } else if ( LQUERY ) {
         RETURN;
      }

      // Quick return if possible

      if ( M == 0 ) {
         RETURN;
      } else if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU( I ) = ZERO;
         } // 10
         RETURN;
      }

      NBMIN = 2;
      NX = 1;
      IWS = M;
      if ( NB > 1 && NB < M ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = MAX( 0, ILAENV( 3, 'ZGERQF', ' ', M, N, -1, -1 ) );
         if ( NX < M ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = M;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               NB = LWORK / LDWORK;
               NBMIN = MAX( 2, ILAENV( 2, 'ZGERQF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < M && NX < M ) {

         // Use blocked code initially.
         // The last kk rows are handled by the block method.

         M1 = MIN( M+1, N );
         KI = ( ( M-NX-1 ) / NB )*NB;
         KK = MIN( M, KI+NB );

         DO 20 I = M - KK + KI + 1, M - KK + 1, -NB;
            IB = MIN( M-I+1, NB );

            // Compute the TZ factorization of the current block
            // A(i:i+ib-1,i:n)

            zlatrz(IB, N-I+1, N-M, A( I, I ), LDA, TAU( I ), WORK );
            if ( I > 1 ) {

               // Form the triangular factor of the block reflector
               // H = H(i+ib-1) . . . H(i+1) H(i)

               zlarzt('Backward', 'Rowwise', N-M, IB, A( I, M1 ), LDA, TAU( I ), WORK, LDWORK );

               // Apply H to A(1:i-1,i:n) from the right

               zlarzb('Right', 'No transpose', 'Backward', 'Rowwise', I-1, N-I+1, IB, N-M, A( I, M1 ), LDA, WORK, LDWORK, A( 1, I ), LDA, WORK( IB+1 ), LDWORK );
            }
         } // 20
         MU = I + NB - 1;
      } else {
         MU = M;
      }

      // Use unblocked code to factor the last or only block

      if (MU > 0) CALL ZLATRZ( MU, N, N-M, A, LDA, TAU, WORK );

      WORK( 1 ) = LWKOPT;

      RETURN;

      // End of ZTZRZF

      }
