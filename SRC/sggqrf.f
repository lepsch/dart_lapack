      SUBROUTINE SGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQRF, SGERQF, SORMQR, XERBLA
*     ..
*     .. External Functions ..
      int                ILAENV
      EXTERNAL           ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           SROUNDUP_LWORK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NB1 = ILAENV( 1, 'SGEQRF', ' ', N, M, -1, -1 )
      NB2 = ILAENV( 1, 'SGERQF', ' ', N, P, -1, -1 )
      NB3 = ILAENV( 1, 'SORMQR', ' ', N, M, P, -1 )
      NB = MAX( NB1, NB2, NB3 )
      LWKOPT = MAX( 1, MAX( N, M, P )*NB )
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, N, M, P ) .AND. .NOT.LQUERY ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     QR factorization of N-by-M matrix A: A = Q*R
*
      CALL SGEQRF( N, M, A, LDA, TAUA, WORK, LWORK, INFO )
      LOPT = INT( WORK( 1 ) )
*
*     Update B := Q**T*B.
*
      CALL SORMQR( 'Left', 'Transpose', N, P, MIN( N, M ), A, LDA, TAUA, B, LDB, WORK, LWORK, INFO )
      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )
*
*     RQ factorization of N-by-P matrix B: B = T*Z.
*
      CALL SGERQF( N, P, B, LDB, TAUB, WORK, LWORK, INFO )
      LWKOPT = MAX( LOPT, INT( WORK( 1 ) ) )
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      RETURN
*
*     End of SGGQRF
*
      END
