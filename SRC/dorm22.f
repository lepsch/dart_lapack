      void dorm22(SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                M, N, N1, N2, LDQ, LDC, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      double             Q( LDQ, * ), C( LDC, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;

      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      int                I, LDWORK, LEN, LWKOPT, NB, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LEFT = LSAME( SIDE, 'L' );
      NOTRAN = LSAME( TRANS, 'N' );
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
      if ( !LEFT && !LSAME( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !LSAME( TRANS, 'N' ) && !LSAME( TRANS, 'T' ) ) {
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
         WORK( 1 ) = DBLE( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('DORM22', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         WORK( 1 ) = 1;
         return;
      }

      // Degenerate cases (N1 = 0 or N2 = 0) are handled using DTRMM.

      if ( N1 == 0 ) {
         dtrmm(SIDE, 'Upper', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK( 1 ) = ONE;
         return;
      } else if ( N2 == 0 ) {
         dtrmm(SIDE, 'Lower', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK( 1 ) = ONE;
         return;
      }

      // Compute the largest chunk size available from the workspace.

      NB = max( 1, min( LWORK, LWKOPT ) / NQ );

      if ( LEFT ) {
         if ( NOTRAN ) {
            DO I = 1, N, NB;
               LEN = min( NB, N-I+1 );
               LDWORK = M;

               // Multiply bottom part of C by Q12.

               dlacpy('All', N1, LEN, C( N2+1, I ), LDC, WORK, LDWORK );
               dtrmm('Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11.

               dgemm('No Transpose', 'No Transpose', N1, LEN, N2, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q21.

               dlacpy('All', N2, LEN, C( 1, I ), LDC, WORK( N1+1 ), LDWORK );
               dtrmm('Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK( N1+1 ), LDWORK );

               // Multiply bottom part of C by Q22.

               dgemm('No Transpose', 'No Transpose', N2, LEN, N1, ONE, Q( N1+1, N2+1 ), LDQ, C( N2+1, I ), LDC, ONE, WORK( N1+1 ), LDWORK );

               // Copy everything back.

               dlacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         } else {
            DO I = 1, N, NB;
               LEN = min( NB, N-I+1 );
               LDWORK = M;

               // Multiply bottom part of C by Q21**T.

               dlacpy('All', N2, LEN, C( N1+1, I ), LDC, WORK, LDWORK );
               dtrmm('Left', 'Upper', 'Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11**T.

               dgemm('Transpose', 'No Transpose', N2, LEN, N1, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q12**T.

               dlacpy('All', N1, LEN, C( 1, I ), LDC, WORK( N2+1 ), LDWORK );
               dtrmm('Left', 'Lower', 'Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK( N2+1 ), LDWORK );

               // Multiply bottom part of C by Q22**T.

               dgemm('Transpose', 'No Transpose', N1, LEN, N2, ONE, Q( N1+1, N2+1 ), LDQ, C( N1+1, I ), LDC, ONE, WORK( N2+1 ), LDWORK );

               // Copy everything back.

               dlacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         }
      } else {
         if ( NOTRAN ) {
            DO I = 1, M, NB;
               LEN = min( NB, M-I+1 );
               LDWORK = LEN;

               // Multiply right part of C by Q21.

               dlacpy('All', LEN, N2, C( I, N1+1 ), LDC, WORK, LDWORK );
               dtrmm('Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11.

               dgemm('No Transpose', 'No Transpose', LEN, N2, N1, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q12.

               dlacpy('All', LEN, N1, C( I, 1 ), LDC, WORK( 1 + N2*LDWORK ), LDWORK );
               dtrmm('Right', 'Lower', 'No Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK( 1 + N2*LDWORK ), LDWORK );

               // Multiply right part of C by Q22.

               dgemm('No Transpose', 'No Transpose', LEN, N1, N2, ONE, C( I, N1+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N2*LDWORK ), LDWORK );

               // Copy everything back.

               dlacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         } else {
            DO I = 1, M, NB;
               LEN = min( NB, M-I+1 );
               LDWORK = LEN;

               // Multiply right part of C by Q12**T.

               dlacpy('All', LEN, N1, C( I, N2+1 ), LDC, WORK, LDWORK );
               dtrmm('Right', 'Lower', 'Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11**T.

               dgemm('No Transpose', 'Transpose', LEN, N1, N2, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q21**T.

               dlacpy('All', LEN, N2, C( I, 1 ), LDC, WORK( 1 + N1*LDWORK ), LDWORK );
               dtrmm('Right', 'Upper', 'Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK( 1 + N1*LDWORK ), LDWORK );

               // Multiply right part of C by Q22**T.

               dgemm('No Transpose', 'Transpose', LEN, N2, N1, ONE, C( I, N2+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N1*LDWORK ), LDWORK );

               // Copy everything back.

               dlacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         }
      }

      WORK( 1 ) = DBLE( LWKOPT );
      return;
      }
