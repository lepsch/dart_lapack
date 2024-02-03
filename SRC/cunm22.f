      SUBROUTINE CUNM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                M, N, N1, N2, LDQ, LDC, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      COMPLEX            Q( LDQ, * ), C( LDC, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;

      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      int                I, LDWORK, LEN, LWKOPT, NB, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      // NQ is the order of Q;
      // NW is the minimum dimension of WORK.

      if ( LEFT ) {
         NQ = M
      } else {
         NQ = N
      }
      NW = NQ
      IF( N1.EQ.0 .OR. N2.EQ.0 ) NW = 1
      if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( N1.LT.0 .OR. N1+N2.NE.NQ ) {
         INFO = -5
      } else if ( N2.LT.0 ) {
         INFO = -6
      } else if ( LDQ.LT.MAX( 1, NQ ) ) {
         INFO = -8
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LWORK.LT.NW .AND. .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO.EQ.0 ) {
         LWKOPT = M*N
         WORK( 1 ) = CMPLX( LWKOPT )
      }

      if ( INFO.NE.0 ) {
         xerbla('CUNM22', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      // Degenerate cases (N1 = 0 or N2 = 0) are handled using CTRMM.

      if ( N1.EQ.0 ) {
         ctrmm(SIDE, 'Upper', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK( 1 ) = ONE
         RETURN
      } else if ( N2.EQ.0 ) {
         ctrmm(SIDE, 'Lower', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC );
         WORK( 1 ) = ONE
         RETURN
      }

      // Compute the largest chunk size available from the workspace.

      NB = MAX( 1, MIN( LWORK, LWKOPT ) / NQ )

      if ( LEFT ) {
         if ( NOTRAN ) {
            DO I = 1, N, NB
               LEN = MIN( NB, N-I+1 )
               LDWORK = M

               // Multiply bottom part of C by Q12.

               clacpy('All', N1, LEN, C( N2+1, I ), LDC, WORK, LDWORK )                CALL CTRMM( 'Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11.

               cgemm('No Transpose', 'No Transpose', N1, LEN, N2, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q21.

               clacpy('All', N2, LEN, C( 1, I ), LDC, WORK( N1+1 ), LDWORK )                CALL CTRMM( 'Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK( N1+1 ), LDWORK );

               // Multiply bottom part of C by Q22.

               cgemm('No Transpose', 'No Transpose', N2, LEN, N1, ONE, Q( N1+1, N2+1 ), LDQ, C( N2+1, I ), LDC, ONE, WORK( N1+1 ), LDWORK );

               // Copy everything back.

               clacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         } else {
            DO I = 1, N, NB
               LEN = MIN( NB, N-I+1 )
               LDWORK = M

               // Multiply bottom part of C by Q21**H.

               clacpy('All', N2, LEN, C( N1+1, I ), LDC, WORK, LDWORK )                CALL CTRMM( 'Left', 'Upper', 'Conjugate', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply top part of C by Q11**H.

               cgemm('Conjugate', 'No Transpose', N2, LEN, N1, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK );

               // Multiply top part of C by Q12**H.

               clacpy('All', N1, LEN, C( 1, I ), LDC, WORK( N2+1 ), LDWORK )                CALL CTRMM( 'Left', 'Lower', 'Conjugate', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK( N2+1 ), LDWORK );

               // Multiply bottom part of C by Q22**H.

               cgemm('Conjugate', 'No Transpose', N1, LEN, N2, ONE, Q( N1+1, N2+1 ), LDQ, C( N1+1, I ), LDC, ONE, WORK( N2+1 ), LDWORK );

               // Copy everything back.

               clacpy('All', M, LEN, WORK, LDWORK, C( 1, I ), LDC );
            }
         }
      } else {
         if ( NOTRAN ) {
            DO I = 1, M, NB
               LEN = MIN( NB, M-I+1 )
               LDWORK = LEN

               // Multiply right part of C by Q21.

               clacpy('All', LEN, N2, C( I, N1+1 ), LDC, WORK, LDWORK )                CALL CTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11.

               cgemm('No Transpose', 'No Transpose', LEN, N2, N1, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q12.

               clacpy('All', LEN, N1, C( I, 1 ), LDC, WORK( 1 + N2*LDWORK ), LDWORK )                CALL CTRMM( 'Right', 'Lower', 'No Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK( 1 + N2*LDWORK ), LDWORK );

               // Multiply right part of C by Q22.

               cgemm('No Transpose', 'No Transpose', LEN, N1, N2, ONE, C( I, N1+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N2*LDWORK ), LDWORK );

               // Copy everything back.

               clacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         } else {
            DO I = 1, M, NB
               LEN = MIN( NB, M-I+1 )
               LDWORK = LEN

               // Multiply right part of C by Q12**H.

               clacpy('All', LEN, N1, C( I, N2+1 ), LDC, WORK, LDWORK )                CALL CTRMM( 'Right', 'Lower', 'Conjugate', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK );

               // Multiply left part of C by Q11**H.

               cgemm('No Transpose', 'Conjugate', LEN, N1, N2, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK );

               // Multiply left part of C by Q21**H.

               clacpy('All', LEN, N2, C( I, 1 ), LDC, WORK( 1 + N1*LDWORK ), LDWORK )                CALL CTRMM( 'Right', 'Upper', 'Conjugate', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK( 1 + N1*LDWORK ), LDWORK );

               // Multiply right part of C by Q22**H.

               cgemm('No Transpose', 'Conjugate', LEN, N2, N1, ONE, C( I, N2+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N1*LDWORK ), LDWORK );

               // Copy everything back.

               clacpy('All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC );
            }
         }
      }

      WORK( 1 ) = CMPLX( LWKOPT )
      RETURN

      // End of CUNM22

      }
