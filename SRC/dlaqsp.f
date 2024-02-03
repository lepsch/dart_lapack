      SUBROUTINE DLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      int                N
      DOUBLE PRECISION   AMAX, SCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J, JC
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
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
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
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
*           Upper triangle of A is stored.
*
            JC = 1
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J
                  AP( JC+I-1 ) = CJ*S( I )*AP( JC+I-1 )
   10          CONTINUE
               JC = JC + J
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            JC = 1
            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, N
                  AP( JC+I-J ) = CJ*S( I )*AP( JC+I-J )
   30          CONTINUE
               JC = JC + N - J + 1
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of DLAQSP
*
      END
