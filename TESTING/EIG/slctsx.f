      LOGICAL          FUNCTION SLCTSX( AR, AI, BETA )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               AI, AR, BETA
*     ..
*
*  =====================================================================
*
*     .. Scalars in Common ..
      LOGICAL            FS
      INTEGER            I, M, MPLUSN, N
*     ..
*     .. Common blocks ..
      COMMON             / MN / M, N, MPLUSN, I, FS
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      IF( FS ) THEN
         I = I + 1
         IF( I.LE.M ) THEN
            SLCTSX = .FALSE.
         ELSE
            SLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            SLCTSX = .TRUE.
         ELSE
            SLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF
*
*       IF( AR/BETA.GT.0.0 )THEN
*          SLCTSX = .TRUE.
*       ELSE
*          SLCTSX = .FALSE.
*       END IF
*
      RETURN
*
*     End of SLCTSX
*
      END