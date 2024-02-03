      SUBROUTINE SPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      REAL               AP( * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                I, JJ
      REAL               SMIN
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPPEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Initialize SMIN and AMAX.
*
      S( 1 ) = AP( 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )
*
      IF( UPPER ) THEN
*
*        UPLO = 'U':  Upper triangle of A is stored.
*        Find the minimum and maximum diagonal elements.
*
         JJ = 1
         DO 10 I = 2, N
            JJ = JJ + I
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   10    CONTINUE
*
      ELSE
*
*        UPLO = 'L':  Lower triangle of A is stored.
*        Find the minimum and maximum diagonal elements.
*
         JJ = 1
         DO 20 I = 2, N
            JJ = JJ + N - I + 2
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   20    CONTINUE
      END IF
*
      IF( SMIN.LE.ZERO ) THEN
*
*        Find the first non-positive diagonal element and return.
*
         DO 30 I = 1, N
            IF( S( I ).LE.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   30    CONTINUE
      ELSE
*
*        Set the scale factors to the reciprocals
*        of the diagonal elements.
*
         DO 40 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   40    CONTINUE
*
*        Compute SCOND = min(S(I)) / max(S(I))
*
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF
      RETURN
*
*     End of SPPEQU
*
      END
