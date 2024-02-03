      SUBROUTINE DSTECH( N, A, B, EIG, TOL, WORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             TOL;
      // ..
      // .. Array Arguments ..
      double             A( * ), B( * ), EIG( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT;
      double             EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSTECT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Check input parameters

      INFO = 0
      IF( N.EQ.0 ) RETURN
      if ( N.LT.0 ) {
         INFO = -1
         RETURN
      }
      if ( TOL.LT.ZERO ) {
         INFO = -5
         RETURN
      }

      // Get machine constants

      EPS = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      UNFLEP = DLAMCH( 'Safe minimum' ) / EPS
      EPS = TOL*EPS

      // Compute maximum absolute eigenvalue, error tolerance

      MX = ABS( EIG( 1 ) )
      for (I = 2; I <= N; I++) { // 10
         MX = MAX( MX, ABS( EIG( I ) ) )
   10 CONTINUE
      EPS = MAX( EPS*MX, UNFLEP )

      // Sort eigenvalues from EIG into WORK

      for (I = 1; I <= N; I++) { // 20
         WORK( I ) = EIG( I )
   20 CONTINUE
      DO 40 I = 1, N - 1
         ISUB = 1
         EMIN = WORK( 1 )
         DO 30 J = 2, N + 1 - I
            if ( WORK( J ).LT.EMIN ) {
               ISUB = J
               EMIN = WORK( J )
            }
   30    CONTINUE
         if ( ISUB.NE.N+1-I ) {
            WORK( ISUB ) = WORK( N+1-I )
            WORK( N+1-I ) = EMIN
         }
   40 CONTINUE

      // TPNT points to singular value at right endpoint of interval
      // BPNT points to singular value at left  endpoint of interval

      TPNT = 1
      BPNT = 1

      // Begin loop over all intervals

   50 CONTINUE
      UPPER = WORK( TPNT ) + EPS
      LOWER = WORK( BPNT ) - EPS

      // Begin loop merging overlapping intervals

   60 CONTINUE
      IF( BPNT.EQ.N ) GO TO 70
      TUPPR = WORK( BPNT+1 ) + EPS
      IF( TUPPR.LT.LOWER ) GO TO 70

      // Merge

      BPNT = BPNT + 1
      LOWER = WORK( BPNT ) - EPS
      GO TO 60
   70 CONTINUE

      // Count singular values in interval [ LOWER, UPPER ]

      dstect(N, A, B, LOWER, NUML );
      dstect(N, A, B, UPPER, NUMU );
      COUNT = NUMU - NUML
      if ( COUNT.NE.BPNT-TPNT+1 ) {

         // Wrong number of singular values in interval

         INFO = TPNT
         GO TO 80
      }
      TPNT = BPNT + 1
      BPNT = TPNT
      IF( TPNT.LE.N ) GO TO 50
   80 CONTINUE
      RETURN

      // End of DSTECH

      }
