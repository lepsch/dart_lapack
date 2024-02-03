      SUBROUTINE ZUNM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                M, N, N1, N2, LDQ, LDC, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         Q( LDQ, * ), C( LDC, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )

      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      int                I, LDWORK, LEN, LWKOPT, NB, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      // NQ is the order of Q;
      // NW is the minimum dimension of WORK.

      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      NW = NQ
      IF( N1.EQ.0 .OR. N2.EQ.0 ) NW = 1
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( N1.LT.0 .OR. N1+N2.NE.NQ ) THEN
         INFO = -5
      ELSE IF( N2.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.MAX( 1, NQ ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF

      IF( INFO.EQ.0 ) THEN
         LWKOPT = M*N
         WORK( 1 ) = DCMPLX( LWKOPT )
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM22', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF

      // Degenerate cases (N1 = 0 or N2 = 0) are handled using ZTRMM.

      IF( N1.EQ.0 ) THEN
         CALL ZTRMM( SIDE, 'Upper', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC )
         WORK( 1 ) = ONE
         RETURN
      ELSE IF( N2.EQ.0 ) THEN
         CALL ZTRMM( SIDE, 'Lower', TRANS, 'Non-Unit', M, N, ONE, Q, LDQ, C, LDC )
         WORK( 1 ) = ONE
         RETURN
      END IF

      // Compute the largest chunk size available from the workspace.

      NB = MAX( 1, MIN( LWORK, LWKOPT ) / NQ )

      IF( LEFT ) THEN
         IF( NOTRAN ) THEN
            DO I = 1, N, NB
               LEN = MIN( NB, N-I+1 )
               LDWORK = M

               // Multiply bottom part of C by Q12.

               CALL ZLACPY( 'All', N1, LEN, C( N2+1, I ), LDC, WORK, LDWORK )                CALL ZTRMM( 'Left', 'Lower', 'No Transpose', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK )

               // Multiply top part of C by Q11.

               CALL ZGEMM( 'No Transpose', 'No Transpose', N1, LEN, N2, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK )

               // Multiply top part of C by Q21.

               CALL ZLACPY( 'All', N2, LEN, C( 1, I ), LDC, WORK( N1+1 ), LDWORK )                CALL ZTRMM( 'Left', 'Upper', 'No Transpose', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK( N1+1 ), LDWORK )

               // Multiply bottom part of C by Q22.

               CALL ZGEMM( 'No Transpose', 'No Transpose', N2, LEN, N1, ONE, Q( N1+1, N2+1 ), LDQ, C( N2+1, I ), LDC, ONE, WORK( N1+1 ), LDWORK )

               // Copy everything back.

               CALL ZLACPY( 'All', M, LEN, WORK, LDWORK, C( 1, I ), LDC )
            END DO
         ELSE
            DO I = 1, N, NB
               LEN = MIN( NB, N-I+1 )
               LDWORK = M

               // Multiply bottom part of C by Q21**H.

               CALL ZLACPY( 'All', N2, LEN, C( N1+1, I ), LDC, WORK, LDWORK )                CALL ZTRMM( 'Left', 'Upper', 'Conjugate', 'Non-Unit', N2, LEN, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK )

               // Multiply top part of C by Q11**H.

               CALL ZGEMM( 'Conjugate', 'No Transpose', N2, LEN, N1, ONE, Q, LDQ, C( 1, I ), LDC, ONE, WORK, LDWORK )

               // Multiply top part of C by Q12**H.

               CALL ZLACPY( 'All', N1, LEN, C( 1, I ), LDC, WORK( N2+1 ), LDWORK )                CALL ZTRMM( 'Left', 'Lower', 'Conjugate', 'Non-Unit', N1, LEN, ONE, Q( 1, N2+1 ), LDQ, WORK( N2+1 ), LDWORK )

               // Multiply bottom part of C by Q22**H.

               CALL ZGEMM( 'Conjugate', 'No Transpose', N1, LEN, N2, ONE, Q( N1+1, N2+1 ), LDQ, C( N1+1, I ), LDC, ONE, WORK( N2+1 ), LDWORK )

               // Copy everything back.

               CALL ZLACPY( 'All', M, LEN, WORK, LDWORK, C( 1, I ), LDC )
            END DO
         END IF
      ELSE
         IF( NOTRAN ) THEN
            DO I = 1, M, NB
               LEN = MIN( NB, M-I+1 )
               LDWORK = LEN

               // Multiply right part of C by Q21.

               CALL ZLACPY( 'All', LEN, N2, C( I, N1+1 ), LDC, WORK, LDWORK )                CALL ZTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK, LDWORK )

               // Multiply left part of C by Q11.

               CALL ZGEMM( 'No Transpose', 'No Transpose', LEN, N2, N1, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK )

               // Multiply left part of C by Q12.

               CALL ZLACPY( 'All', LEN, N1, C( I, 1 ), LDC, WORK( 1 + N2*LDWORK ), LDWORK )                CALL ZTRMM( 'Right', 'Lower', 'No Transpose', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK( 1 + N2*LDWORK ), LDWORK )

               // Multiply right part of C by Q22.

               CALL ZGEMM( 'No Transpose', 'No Transpose', LEN, N1, N2, ONE, C( I, N1+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N2*LDWORK ), LDWORK )

               // Copy everything back.

               CALL ZLACPY( 'All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC )
            END DO
         ELSE
            DO I = 1, M, NB
               LEN = MIN( NB, M-I+1 )
               LDWORK = LEN

               // Multiply right part of C by Q12**H.

               CALL ZLACPY( 'All', LEN, N1, C( I, N2+1 ), LDC, WORK, LDWORK )                CALL ZTRMM( 'Right', 'Lower', 'Conjugate', 'Non-Unit', LEN, N1, ONE, Q( 1, N2+1 ), LDQ, WORK, LDWORK )

               // Multiply left part of C by Q11**H.

               CALL ZGEMM( 'No Transpose', 'Conjugate', LEN, N1, N2, ONE, C( I, 1 ), LDC, Q, LDQ, ONE, WORK, LDWORK )

               // Multiply left part of C by Q21**H.

               CALL ZLACPY( 'All', LEN, N2, C( I, 1 ), LDC, WORK( 1 + N1*LDWORK ), LDWORK )                CALL ZTRMM( 'Right', 'Upper', 'Conjugate', 'Non-Unit', LEN, N2, ONE, Q( N1+1, 1 ), LDQ, WORK( 1 + N1*LDWORK ), LDWORK )

               // Multiply right part of C by Q22**H.

               CALL ZGEMM( 'No Transpose', 'Conjugate', LEN, N2, N1, ONE, C( I, N2+1 ), LDC, Q( N1+1, N2+1 ), LDQ, ONE, WORK( 1 + N1*LDWORK ), LDWORK )

               // Copy everything back.

               CALL ZLACPY( 'All', LEN, N, WORK, LDWORK, C( I, 1 ), LDC )
            END DO
         END IF
      END IF

      WORK( 1 ) = DCMPLX( LWKOPT )
      RETURN

      // End of ZUNM22

      END
