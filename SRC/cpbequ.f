      SUBROUTINE CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               S( * )
      COMPLEX            AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J;
      REAL               SMIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPBEQU', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF

      IF( UPPER ) THEN
         J = KD + 1
      } else {
         J = 1
      END IF

      // Initialize SMIN and AMAX.

      S( 1 ) = REAL( AB( J, 1 ) )
      SMIN = S( 1 )
      AMAX = S( 1 )

      // Find the minimum and maximum diagonal elements.

      DO 10 I = 2, N
         S( I ) = REAL( AB( J, I ) )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE

      IF( SMIN.LE.ZERO ) THEN

         // Find the first non-positive diagonal element and return.

         DO 20 I = 1, N
            IF( S( I ).LE.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   20    CONTINUE
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         DO 30 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   30    CONTINUE

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF
      RETURN

      // End of CPBEQU

      }
