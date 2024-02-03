      SUBROUTINE DLARRR( N, D, E, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, INFO;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      // ..


*  =====================================================================

      // .. Parameters ..
      double             ZERO, RELCOND;
      const              ZERO = 0.0D0, RELCOND = 0.999D0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      bool               YESREL;
      double             EPS, SAFMIN, SMLNUM, RMIN, TMP, TMP2, OFFDIG, OFFDIG2;

      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         INFO = 0
         RETURN
      }

      // As a default, do NOT go for relative-accuracy preserving computations.
      INFO = 1

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      RMIN = SQRT( SMLNUM )

      // Tests for relative accuracy

      // Test for scaled diagonal dominance
      // Scale the diagonal entries to one and check whether the sum of the
      // off-diagonals is less than one

      // The sdd relative error bounds have a 1/(1- 2*x) factor in them,
      // x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
      // accuracy is promised.  In the notation of the code fragment below,
      // 1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
      // We don't think it is worth going into "sdd mode" unless the relative
      // condition number is reasonable, not 1/macheps.
      // The threshold should be compatible with other thresholds used in the
      // code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
      // to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
      // instead of the current OFFDIG + OFFDIG2 < 1

      YESREL = .TRUE.
      OFFDIG = ZERO
      TMP = SQRT(ABS(D(1)))
      IF (TMP.LT.RMIN) YESREL = .FALSE.
      IF(.NOT.YESREL) GOTO 11
      for (I = 2; I <= N; I++) { // 10
         TMP2 = SQRT(ABS(D(I)))
         IF (TMP2.LT.RMIN) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         OFFDIG2 = ABS(E(I-1))/(TMP*TMP2)
         IF(OFFDIG+OFFDIG2.GE.RELCOND) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         TMP = TMP2
         OFFDIG = OFFDIG2
      } // 10
      } // 11

      if ( YESREL ) {
         INFO = 0
         RETURN
      } else {
      ENDIF



      // *** MORE TO BE IMPLEMENTED ***



      // Test if the lower bidiagonal matrix L from T = L D L^T
      // (zero shift facto) is well conditioned



      // Test if the upper bidiagonal matrix U from T = U D U^T
      // (zero shift facto) is well conditioned.
      // In this case, the matrix needs to be flipped and, at the end
      // of the eigenvector computation, the flip needs to be applied
      // to the computed eigenvectors (and the support)



      RETURN

      // End of DLARRR

      }
