      SUBROUTINE DSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER, LQUERY
      int                LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSYTRI_3X, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     Determine the block size
*
      IF( N.EQ.0 ) THEN
         LWKOPT = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'DSYTRI_3', UPLO, N, -1, -1, -1 ) )
         LWKOPT = ( N+NB+1 ) * ( NB+3 )
      END IF
      WORK( 1 ) = LWKOPT
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.LWKOPT .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRI_3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      CALL DSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYTRI_3
*
      END
