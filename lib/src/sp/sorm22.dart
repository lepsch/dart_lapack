      void sorm22(final int SIDE, final int TRANS, final int M, final int N, final int N1, final int N2, final Matrix<double> Q, final int LDQ, final Matrix<double> C, final int LDC, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      String             SIDE, TRANS;
      int                M, N, N1, N2, LDQ, LDC, LWORK, INFO;
      double               Q( LDQ, * ), C( LDC, * ), WORK( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;

      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      int                I, LDWORK, LEN, LWKOPT, NB, NQ, NW;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, STRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      LEFT = lsame( SIDE, 'L' );
      NOTRAN = lsame( TRANS, 'N' );
      LQUERY = ( LWORK == -1 );

      // NQ is the order of Q;
      // NW is the minimum dimension of WORK.

      if ( LEFT ) {
         NQ = M;
      } else {
         NQ = N;
      }
      NW = NQ;
      if (N1 == 0 || N2 == 0) NW = 1;
      if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'T' ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( N1 < 0 || N1+N2 != NQ ) {
         INFO = -5;
      } else if ( N2 < 0 ) {
         INFO = -6;
      } else if ( LDQ < max( 1, NQ ) ) {
         INFO = -8;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      } else if ( LWORK < NW && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO == 0 ) {
         LWKOPT = M*N;
         WORK[1] = SROUNDUP_LWORK( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('SORM22', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         WORK[1] = 1;
         return;
      }

      // Degenerate cases (N1 = 0 or N2 = 0) are handled using STRMM.

      if ( N1 == 0 ) {
         strmm(SIDE, 'Upper', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK[1] = ONE;
         return;
      } else if ( N2 == 0 ) {
         strmm(SIDE, 'Lower', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK[1] = ONE;
         return;
      }

      // Compute the largest chunk size available from the workspace.

      NB = max( 1, min( LWORK, LWKOPT ) / NQ );

      if ( LEFT ) {
         if ( NOTRAN ) {
            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
               LEN = min( NB, N-I+1 );
               LDWORK = M;

               // Multiply bottom part of C by Q12.

               slacpy('All', N1, LEN, C( N2+1, I ), LDC, WORK, LDWORK );
               strmm('Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11.

               sgemm('No Transpose', 'No Transpose', N1, LEN, N2, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q21.

               slacpy('All', N2, LEN, C( 1, I ), LDC, WORK( N1+1 ), LDWORK );
               strmm('Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK( N1+1 ), LDWORK );

               // Multiply bottom part of C by Q22.

               sgemm('No Transpose', 'No Transpose', N2, LEN, N1, ONE, Q( N1+1, N2+1 ), LDQ, C( N2+1, I ), LDC, ONE, WORK( N1+1 ), LDWORK );

               // Copy everything back.

               slacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         } else {
            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) {
               LEN = min( NB, N-I+1 );
               LDWORK = M;

               // Multiply bottom part of C by Q21**T.

               slacpy('All', N2, LEN, C( N1+1, I ), LDC, WORK, LDWORK );
               strmm('Left', 'Upper', 'Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11**T.

               sgemm('Transpose', 'No Transpose', N2, LEN, N1, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q12**T.

               slacpy('All', N1, LEN, C( 1, I ), LDC, WORK( N2+1 ), LDWORK );
               strmm('Left', 'Lower', 'Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK( N2+1 ), LDWORK );

               // Multiply bottom part of C by Q22**T.

               sgemm('Transpose', 'No Transpose', N1, LEN, N2, ONE, Q( N1+1, N2+1 ), LDQ, C( N1+1, I ), LDC, ONE, WORK( N2+1 ), LDWORK );

               // Copy everything back.

               slacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         }
      } else {
         if ( NOTRAN ) {
            for (I = 1; NB < 0 ? I >= M : I <= M; I += NB) {
               LEN = min( NB, M-I+1 );
               LDWORK = LEN;

               // Multiply right part of C by Q21.

               slacpy('All', LEN, N2, C( I, N1+1 ), LDC, WORK, LDWORK );
               strmm('Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11.

               sgemm('No Transpose', 'No Transpose', LEN, N2, N1, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q12.

               slacpy('All', LEN, N1, C( I, 1 ), LDC, WORK( 1 + N2*LDWORK ), LDWORK );
               strmm('Right', 'Lower', 'No Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK( 1 + N2*LDWORK ), LDWORK );

               // Multiply right part of C by Q22.

               sgemm('No Transpose', 'No Transpose', LEN, N1, N2, ONE, C( I, N1+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N2*LDWORK ), LDWORK );

               // Copy everything back.

               slacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         } else {
            for (I = 1; NB < 0 ? I >= M : I <= M; I += NB) {
               LEN = min( NB, M-I+1 );
               LDWORK = LEN;

               // Multiply right part of C by Q12**T.

               slacpy('All', LEN, N1, C( I, N2+1 ), LDC, WORK, LDWORK );
               strmm('Right', 'Lower', 'Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11**T.

               sgemm('No Transpose', 'Transpose', LEN, N1, N2, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q21**T.

               slacpy('All', LEN, N2, C( I, 1 ), LDC, WORK( 1 + N1*LDWORK ), LDWORK );
               strmm('Right', 'Upper', 'Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK( 1 + N1*LDWORK ), LDWORK );

               // Multiply right part of C by Q22**T.

               sgemm('No Transpose', 'Transpose', LEN, N2, N1, ONE, C( I, N2+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N1*LDWORK ), LDWORK );

               // Copy everything back.

               slacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         }
      }

      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      }
