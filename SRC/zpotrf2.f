      RECURSIVE SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = (1.0D+0, 0.0D+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            N1, N2, IINFO
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHERK, ZTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTRF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     N=1 case
*
      IF( N.EQ.1 ) THEN
*
*        Test for non-positive-definiteness
*
         AJJ = DBLE( A( 1, 1 ) )
         IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Factor
*
         A( 1, 1 ) = SQRT( AJJ )
*
*     Use recursive code
*
      ELSE
         N1 = N/2
         N2 = N-N1
*
*        Factor A11
*
         CALL ZPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF
*
*        Compute the Cholesky factorization A = U**H*U
*
         IF( UPPER ) THEN
*
*           Update and scale A12
*
            CALL ZTRSM( 'L', 'U', 'C', 'N', N1, N2, CONE,
     $                  A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )
*
*           Update and factor A22
*
            CALL ZHERK( UPLO, 'C', N2, N1, -ONE, A( 1, N1+1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL ZPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
*
*        Compute the Cholesky factorization A = L*L**H
*
         ELSE
*
*           Update and scale A21
*
            CALL ZTRSM( 'R', 'L', 'C', 'N', N2, N1, CONE,
     $                  A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
*
*           Update and factor A22
*
            CALL ZHERK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL ZPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZPOTRF2
*
      END