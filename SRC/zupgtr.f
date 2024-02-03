      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDQ, N;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, IJ, J;
*     ..
*     .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA, ZUNG2L, ZUNG2R
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPGTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to ZHPTRD with UPLO = 'U'
*
*        Unpack the vectors which define the elementary reflectors and
*        set the last row and column of Q equal to those of the unit
*        matrix
*
         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = CZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = CZERO
   30    CONTINUE
         Q( N, N ) = CONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL ZUNG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to ZHPTRD with UPLO = 'L'.
*
*        Unpack the vectors which define the elementary reflectors and
*        set the first row and column of Q equal to those of the unit
*        matrix
*
         Q( 1, 1 ) = CONE
         DO 40 I = 2, N
            Q( I, 1 ) = CZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = CZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL ZUNG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, IINFO )
         END IF
      END IF
      RETURN
*
*     End of ZUPGTR
*
      END
