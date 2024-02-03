      int              FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      int                ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN, REAL
*     ..
*     .. External Functions ..
      int                IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Arrays in Common ..
      int                IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Save statement ..
      SAVE               / CLAENV /
*     ..
*     .. Executable Statements ..
*
      IF( ISPEC.GE.1 .AND. ISPEC.LE.5 ) THEN
*
*        Return a value from the common block.
*
         IF ( NAME(2:6).EQ.'GEQR ' ) THEN
            IF (N3.EQ.2) THEN
               ILAENV = IPARMS ( 2 )
            ELSE
               ILAENV = IPARMS ( 1 )
            END IF
         ELSE IF ( NAME(2:6).EQ.'GELQ ' ) THEN
            IF (N3.EQ.2) THEN
               ILAENV = IPARMS ( 2 )
            ELSE
               ILAENV = IPARMS ( 1 )
            END IF
         ELSE
            ILAENV = IPARMS( ISPEC )
         END IF
*
      ELSE IF( ISPEC.EQ.6 ) THEN
*
*        Compute SVD crossover point.
*
         ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
*
      ELSE IF( ISPEC.GE.7 .AND. ISPEC.LE.9 ) THEN
*
*        Return a value from the common block.
*
         ILAENV = IPARMS( ISPEC )
*
      ELSE IF( ISPEC.EQ.10 ) THEN
*
*        IEEE NaN arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ILAENV = 1
         IF( ILAENV.EQ.1 ) THEN
            ILAENV = IEEECK( 1, 0.0, 1.0 )
         END IF
*
      ELSE IF( ISPEC.EQ.11 ) THEN
*
*        Infinity arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ILAENV = 1
         IF( ILAENV.EQ.1 ) THEN
            ILAENV = IEEECK( 0, 0.0, 1.0 )
         END IF
*
      ELSE
*
*        Invalid value for ISPEC
*
         ILAENV = -1
      END IF
*
      RETURN
*
*     End of ILAENV
*
      END
      int     FUNCTION ILAENV2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      int                ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local variables ..
      int                IISPEC
*     .. External Functions ..
      int                IPARAM2STAGE
      EXTERNAL           IPARAM2STAGE
*     ..
*     .. Arrays in Common ..
      int                IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Save statement ..
      SAVE               / CLAENV /
*     ..
*     .. Executable Statements ..
*
      IF(( ISPEC.GE.1 ) .AND. (ISPEC.LE.5)) THEN
*
*     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
*
         IF( ISPEC.EQ.1 ) THEN
             ILAENV2STAGE = IPARMS( 1 )
         ELSE
             IISPEC = 16 + ISPEC
             ILAENV2STAGE = IPARAM2STAGE( IISPEC, NAME, OPTS, N1, N2, N3, N4 )
         ENDIF
*
      ELSE
*
*        Invalid value for ISPEC
*
         ILAENV2STAGE = -1
      END IF
*
      RETURN
      END
