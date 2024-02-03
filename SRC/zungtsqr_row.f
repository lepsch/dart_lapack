      SUBROUTINE ZUNGTSQR_ROW( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDT, LWORK, M, N, MB, NB
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CONE, CZERO
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            NBLOCAL, MB2, M_PLUS_ONE, ITMP, IB_BOTTOM, LWORKOPT, NUM_ALL_ROW_BLOCKS, JB_T, IB, IMB, KB, KB_LAST, KNB, MB1
*     ..
*     .. Local Arrays ..
      COMPLEX*16         DUMMY( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLARFB_GETT, ZLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      LQUERY  = LWORK.EQ.-1
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.LT.N ) THEN
         INFO = -2
      ELSE IF( MB.LE.N ) THEN
         INFO = -3
      ELSE IF( NB.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDT.LT.MAX( 1, MIN( NB, N ) ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      NBLOCAL = MIN( NB, N )
*
*     Determine the workspace size.
*
      IF( INFO.EQ.0 ) THEN
         LWORKOPT = NBLOCAL * MAX( NBLOCAL, ( N - NBLOCAL ) )
      END IF
*
*     Handle error in the input parameters and handle the workspace query.
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGTSQR_ROW', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
         WORK( 1 ) = DCMPLX( LWORKOPT )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) THEN
         WORK( 1 ) = DCMPLX( LWORKOPT )
         RETURN
      END IF
*
*     (0) Set the upper-triangular part of the matrix A to zero and
*     its diagonal elements to one.
*
      CALL ZLASET('U', M, N, CZERO, CONE, A, LDA )
*
*     KB_LAST is the column index of the last column block reflector
*     in the matrices T and V.
*
      KB_LAST = ( ( N-1 ) / NBLOCAL ) * NBLOCAL + 1
*
*
*     (1) Bottom-up loop over row blocks of A, except the top row block.
*     NOTE: If MB>=M, then the loop is never executed.
*
      IF ( MB.LT.M ) THEN
*
*        MB2 is the row blocking size for the row blocks before the
*        first top row block in the matrix A. IB is the row index for
*        the row blocks in the matrix A before the first top row block.
*        IB_BOTTOM is the row index for the last bottom row block
*        in the matrix A. JB_T is the column index of the corresponding
*        column block in the matrix T.
*
*        Initialize variables.
*
*        NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
*        including the first row block.
*
         MB2 = MB - N
         M_PLUS_ONE = M + 1
         ITMP = ( M - MB - 1 ) / MB2
         IB_BOTTOM = ITMP * MB2 + MB + 1
         NUM_ALL_ROW_BLOCKS = ITMP + 2
         JB_T = NUM_ALL_ROW_BLOCKS * N + 1
*
         DO IB = IB_BOTTOM, MB+1, -MB2
*
*           Determine the block size IMB for the current row block
*           in the matrix A.
*
            IMB = MIN( M_PLUS_ONE - IB, MB2 )
*
*           Determine the column index JB_T for the current column block
*           in the matrix T.
*
            JB_T = JB_T - N
*
*           Apply column blocks of H in the row block from right to left.
*
*           KB is the column index of the current column block reflector
*           in the matrices T and V.
*
            DO KB = KB_LAST, 1, -NBLOCAL
*
*              Determine the size of the current column block KNB in
*              the matrices T and V.
*
               KNB = MIN( NBLOCAL, N - KB + 1 )
*
               CALL ZLARFB_GETT( 'I', IMB, N-KB+1, KNB, T( 1, JB_T+KB-1 ), LDT, A( KB, KB ), LDA, A( IB, KB ), LDA, WORK, KNB )
*
            END DO
*
         END DO
*
      END IF
*
*     (2) Top row block of A.
*     NOTE: If MB>=M, then we have only one row block of A of size M
*     and we work on the entire matrix A.
*
      MB1 = MIN( MB, M )
*
*     Apply column blocks of H in the top row block from right to left.
*
*     KB is the column index of the current block reflector in
*     the matrices T and V.
*
      DO KB = KB_LAST, 1, -NBLOCAL
*
*        Determine the size of the current column block KNB in
*        the matrices T and V.
*
         KNB = MIN( NBLOCAL, N - KB + 1 )
*
         IF( MB1-KB-KNB+1.EQ.0 ) THEN
*
*           In SLARFB_GETT parameters, when M=0, then the matrix B
*           does not exist, hence we need to pass a dummy array
*           reference DUMMY(1,1) to B with LDDUMMY=1.
*
            CALL ZLARFB_GETT( 'N', 0, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, DUMMY( 1, 1 ), 1, WORK, KNB )
         ELSE
            CALL ZLARFB_GETT( 'N', MB1-KB-KNB+1, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, A( KB+KNB, KB), LDA, WORK, KNB )

         END IF
*
      END DO
*
      WORK( 1 ) = DCMPLX( LWORKOPT )
      RETURN
*
*     End of ZUNGTSQR_ROW
*
      END
