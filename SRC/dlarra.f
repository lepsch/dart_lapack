      SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM, NSPLIT, ISPLIT, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N, NSPLIT;
      double              SPLTOL, TNRM;
      // ..
      // .. Array Arguments ..
      int                ISPLIT( * );
      double             D( * ), E( * ), E2( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             EABS, TMP1;

      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      INFO = 0
      NSPLIT = 1

      // Quick return if possible

      IF( N.LE.0 ) THEN
         RETURN
      END IF

      // Compute splitting points
      IF(SPLTOL.LT.ZERO) THEN
         // Criterion based on absolute off-diagonal value
         TMP1 = ABS(SPLTOL)* TNRM
         DO 9 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. TMP1) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I
               NSPLIT = NSPLIT + 1
            END IF
 9       CONTINUE
      } else {
         // Criterion that guarantees relative accuracy
         DO 10 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. SPLTOL * SQRT(ABS(D(I)))*SQRT(ABS(D(I+1))) ) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I
               NSPLIT = NSPLIT + 1
            END IF
 10      CONTINUE
      ENDIF
      ISPLIT( NSPLIT ) = N

      RETURN

      // End of DLARRA

      }
