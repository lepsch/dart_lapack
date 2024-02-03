      SUBROUTINE ZLATSP( UPLO, N, X, ISEED )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N
*     ..
*     .. Array Arguments ..
      int                ISEED( * )
      COMPLEX*16         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         EYE
      PARAMETER          ( EYE = ( 0.0D0, 1.0D0 ) )
*     ..
*     .. Local Scalars ..
      int                J, JJ, N5
      DOUBLE PRECISION   ALPHA, ALPHA3, BETA
      COMPLEX*16         A, B, C, R
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLARND
      EXTERNAL           ZLARND
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Initialize constants
*
      ALPHA = ( 1.D0+SQRT( 17.D0 ) ) / 8.D0
      BETA = ALPHA - 1.D0 / 1000.D0
      ALPHA3 = ALPHA*ALPHA*ALPHA
*
*     Fill the matrix with zeros.
*
      DO 10 J = 1, N*( N+1 ) / 2
         X( J ) = 0.0D0
   10 CONTINUE
*
*     UPLO = 'U':  Upper triangular storage
*
      IF( UPLO.EQ.'U' ) THEN
         N5 = N / 5
         N5 = N - 5*N5 + 1
*
         JJ = N*( N+1 ) / 2
         DO 20 J = N, N5, -5
            A = ALPHA3*ZLARND( 5, ISEED )
            B = ZLARND( 5, ISEED ) / ALPHA
            C = A - 2.D0*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            X( JJ ) = ZLARND( 2, ISEED )
            JJ = JJ - ( J-3 )
            X( JJ ) = ZLARND( 2, ISEED )
            IF( ABS( X( JJ+( J-3 ) ) ).GT.ABS( X( JJ ) ) ) THEN
               X( JJ+( J-4 ) ) = 2.0D0*X( JJ+( J-3 ) )
            ELSE
               X( JJ+( J-4 ) ) = 2.0D0*X( JJ )
            END IF
            JJ = JJ - ( J-4 )
   20    CONTINUE
*
*        Clean-up for N not a multiple of 5.
*
         J = N5 - 1
         IF( J.GT.2 ) THEN
            A = ALPHA3*ZLARND( 5, ISEED )
            B = ZLARND( 5, ISEED ) / ALPHA
            C = A - 2.D0*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            J = J - 3
         END IF
         IF( J.GT.1 ) THEN
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ-J ) = ZLARND( 2, ISEED )
            IF( ABS( X( JJ ) ).GT.ABS( X( JJ-J ) ) ) THEN
               X( JJ-1 ) = 2.0D0*X( JJ )
            ELSE
               X( JJ-1 ) = 2.0D0*X( JJ-J )
            END IF
            JJ = JJ - J - ( J-1 )
            J = J - 2
         ELSE IF( J.EQ.1 ) THEN
            X( JJ ) = ZLARND( 2, ISEED )
            J = J - 1
         END IF
*
*     UPLO = 'L':  Lower triangular storage
*
      ELSE
         N5 = N / 5
         N5 = N5*5
*
         JJ = 1
         DO 30 J = 1, N5, 5
            A = ALPHA3*ZLARND( 5, ISEED )
            B = ZLARND( 5, ISEED ) / ALPHA
            C = A - 2.D0*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            X( JJ ) = ZLARND( 2, ISEED )
            JJ = JJ + ( N-J-2 )
            X( JJ ) = ZLARND( 2, ISEED )
            IF( ABS( X( JJ-( N-J-2 ) ) ).GT.ABS( X( JJ ) ) ) THEN
               X( JJ-( N-J-2 )+1 ) = 2.0D0*X( JJ-( N-J-2 ) )
            ELSE
               X( JJ-( N-J-2 )+1 ) = 2.0D0*X( JJ )
            END IF
            JJ = JJ + ( N-J-3 )
   30    CONTINUE
*
*        Clean-up for N not a multiple of 5.
*
         J = N5 + 1
         IF( J.LT.N-1 ) THEN
            A = ALPHA3*ZLARND( 5, ISEED )
            B = ZLARND( 5, ISEED ) / ALPHA
            C = A - 2.D0*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            J = J + 3
         END IF
         IF( J.LT.N ) THEN
            X( JJ ) = ZLARND( 2, ISEED )
            X( JJ+( N-J+1 ) ) = ZLARND( 2, ISEED )
            IF( ABS( X( JJ ) ).GT.ABS( X( JJ+( N-J+1 ) ) ) ) THEN
               X( JJ+1 ) = 2.0D0*X( JJ )
            ELSE
               X( JJ+1 ) = 2.0D0*X( JJ+( N-J+1 ) )
            END IF
            JJ = JJ + ( N-J+1 ) + ( N-J )
            J = J + 2
         ELSE IF( J.EQ.N ) THEN
            X( JJ ) = ZLARND( 2, ISEED )
            JJ = JJ + ( N-J+1 )
            J = J + 1
         END IF
      END IF
*
      RETURN
*
*     End of ZLATSP
*
      END
