      SUBROUTINE CTPTTR( UPLO, N, AP, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      int                I, J, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPTTR', -INFO )
         RETURN
      END IF
*
      IF( LOWER ) THEN
         K = 0
         DO J = 1, N
            DO I = J, N
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      ELSE
         K = 0
         DO J = 1, N
            DO I = 1, J
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      END IF
*
*
      RETURN
*
*     End of CTPTTR
*
      END
