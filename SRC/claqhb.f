      SUBROUTINE CLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            KD, LDAB, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      REAL               S( * )
      COMPLEX            AB( LDAB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, THRESH
      PARAMETER          ( ONE = 1.0E+0, THRESH = 0.1E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
      ELSE
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A is stored in band format.
*
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = MAX( 1, J-KD ), J - 1
                  AB( KD+1+I-J, J ) = CJ*S( I )*AB( KD+1+I-J, J )
   10          CONTINUE
               AB( KD+1, J ) = CJ*CJ*REAL( AB( KD+1, J ) )
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            DO 40 J = 1, N
               CJ = S( J )
               AB( 1, J ) = CJ*CJ*REAL( AB( 1, J ) )
               DO 30 I = J + 1, MIN( N, J+KD )
                  AB( 1+I-J, J ) = CJ*S( I )*AB( 1+I-J, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of CLAQHB
*
      END