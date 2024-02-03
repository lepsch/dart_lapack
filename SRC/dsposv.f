      SUBROUTINE DSPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, SWORK, ITER, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITER, LDA, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               SWORK( * )
      double             A( LDA, * ), B( LDB, * ), WORK( N, * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      bool               DOITREF;
      const              DOITREF = .TRUE. ;

      int                ITERMAX;
      const              ITERMAX = 30 ;

      double             BWDMAX;
      const              BWDMAX = 1.0E+00 ;

      double             NEGONE, ONE;
      const              NEGONE = -1.0D+0, ONE = 1.0D+0 ;

      // .. Local Scalars ..
      int                I, IITER, PTSA, PTSX;
      double             ANRM, CTE, EPS, RNRM, XNRM;

      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSYMM, DLACPY, DLAT2S, DLAG2S, SLAG2D, SPOTRF, SPOTRS, DPOTRF, DPOTRS, XERBLA
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DLANSY;
      bool               LSAME;
      // EXTERNAL IDAMAX, DLAMCH, DLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0
      ITER = 0

      // Test the input parameters.

      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('DSPOSV', -INFO );
         RETURN
      }

      // Quick return if (N.EQ.0).

      IF( N.EQ.0 ) RETURN

      // Skip single precision iterative refinement if a priori slower
      // than double precision factorization.

      if ( .NOT.DOITREF ) {
         ITER = -1
         GO TO 40
      }

      // Compute some constants.

      ANRM = DLANSY( 'I', UPLO, N, A, LDA, WORK )
      EPS = DLAMCH( 'Epsilon' )
      CTE = ANRM*EPS*SQRT( DBLE( N ) )*BWDMAX

      // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

      PTSA = 1
      PTSX = PTSA + N*N

      // Convert B from double precision to single precision and store the
      // result in SX.

      dlag2s(N, NRHS, B, LDB, SWORK( PTSX ), N, INFO );

      if ( INFO.NE.0 ) {
         ITER = -2
         GO TO 40
      }

      // Convert A from double precision to single precision and store the
      // result in SA.

      dlat2s(UPLO, N, A, LDA, SWORK( PTSA ), N, INFO );

      if ( INFO.NE.0 ) {
         ITER = -2
         GO TO 40
      }

      // Compute the Cholesky factorization of SA.

      spotrf(UPLO, N, SWORK( PTSA ), N, INFO );

      if ( INFO.NE.0 ) {
         ITER = -3
         GO TO 40
      }

      // Solve the system SA*SX = SB.

      spotrs(UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO );

      // Convert SX back to double precision

      slag2d(N, NRHS, SWORK( PTSX ), N, X, LDX, INFO );

      // Compute R = B - AX (R is WORK).

      dlacpy('All', N, NRHS, B, LDB, WORK, N );

      dsymm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER=0 and return.

      for (I = 1; I <= NRHS; I++) {
         XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
         RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
         IF( RNRM.GT.XNRM*CTE ) GO TO 10
      END DO

      // If we are here, the NRHS normwise backward errors satisfy the
      // stopping criterion. We are good to exit.

      ITER = 0
      RETURN

   10 CONTINUE

      for (IITER = 1; IITER <= ITERMAX; IITER++) { // 30

         // Convert R (in WORK) from double precision to single precision
         // and store the result in SX.

         dlag2s(N, NRHS, WORK, N, SWORK( PTSX ), N, INFO );

         if ( INFO.NE.0 ) {
            ITER = -2
            GO TO 40
         }

         // Solve the system SA*SX = SR.

         spotrs(UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO );

         // Convert SX back to double precision and update the current
         // iterate.

         slag2d(N, NRHS, SWORK( PTSX ), N, WORK, N, INFO );

         for (I = 1; I <= NRHS; I++) {
            daxpy(N, ONE, WORK( 1, I ), 1, X( 1, I ), 1 );
         END DO

         // Compute R = B - AX (R is WORK).

         dlacpy('All', N, NRHS, B, LDB, WORK, N );

         dsymm('L', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

         // Check whether the NRHS normwise backward errors satisfy the
         // stopping criterion. If yes, set ITER=IITER>0 and return.

         for (I = 1; I <= NRHS; I++) {
            XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
            RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
            IF( RNRM.GT.XNRM*CTE ) GO TO 20
         END DO

         // If we are here, the NRHS normwise backward errors satisfy the
         // stopping criterion, we are good to exit.

         ITER = IITER

         RETURN

   20    CONTINUE

   30 CONTINUE

      // If we are at this place of the code, this is because we have
      // performed ITER=ITERMAX iterations and never satisfied the
      // stopping criterion, set up the ITER flag accordingly and follow
      // up on double precision routine.

      ITER = -ITERMAX - 1

   40 CONTINUE

      // Single-precision iterative refinement failed to converge to a
      // satisfactory solution, so we resort to double precision.

      dpotrf(UPLO, N, A, LDA, INFO );

      IF( INFO.NE.0 ) RETURN

      dlacpy('All', N, NRHS, B, LDB, X, LDX );
      dpotrs(UPLO, N, NRHS, A, LDA, X, LDX, INFO );

      RETURN

      // End of DSPOSV

      }
