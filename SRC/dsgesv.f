      SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, SWORK, ITER, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, ITER, LDA, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               SWORK( * )
      double             A( LDA, * ), B( LDB, * ), WORK( N, * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      bool               DOITREF;
      const              DOITREF = true ;

      int                ITERMAX;
      const              ITERMAX = 30 ;

      double             BWDMAX;
      const              BWDMAX = 1.0e+00 ;

      double             NEGONE, ONE;
      const              NEGONE = -1.0, ONE = 1.0 ;

      // .. Local Scalars ..
      int                I, IITER, PTSA, PTSX;
      double             ANRM, CTE, EPS, RNRM, XNRM;

      // .. External Subroutines ..
      // EXTERNAL DAXPY, DGEMM, DLACPY, DLAG2S, DGETRF, DGETRS, SGETRF, SGETRS, SLAG2D, XERBLA
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DLANGE;
      // EXTERNAL IDAMAX, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0
      ITER = 0

      // Test the input parameters.

      if ( N < 0 ) {
         INFO = -1
      } else if ( NRHS < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('DSGESV', -INFO );
         RETURN
      }

      // Quick return if (N == 0).

      if (N == 0) RETURN;

      // Skip single precision iterative refinement if a priori slower
      // than double precision factorization.

      if ( .NOT.DOITREF ) {
         ITER = -1
         GO TO 40
      }

      // Compute some constants.

      ANRM = DLANGE( 'I', N, N, A, LDA, WORK )
      EPS = DLAMCH( 'Epsilon' )
      CTE = ANRM*EPS*SQRT( DBLE( N ) )*BWDMAX

      // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

      PTSA = 1
      PTSX = PTSA + N*N

      // Convert B from double precision to single precision and store the
      // result in SX.

      dlag2s(N, NRHS, B, LDB, SWORK( PTSX ), N, INFO );

      if ( INFO != 0 ) {
         ITER = -2
         GO TO 40
      }

      // Convert A from double precision to single precision and store the
      // result in SA.

      dlag2s(N, N, A, LDA, SWORK( PTSA ), N, INFO );

      if ( INFO != 0 ) {
         ITER = -2
         GO TO 40
      }

      // Compute the LU factorization of SA.

      sgetrf(N, N, SWORK( PTSA ), N, IPIV, INFO );

      if ( INFO != 0 ) {
         ITER = -3
         GO TO 40
      }

      // Solve the system SA*SX = SB.

      sgetrs('No transpose', N, NRHS, SWORK( PTSA ), N, IPIV, SWORK( PTSX ), N, INFO );

      // Convert SX back to double precision

      slag2d(N, NRHS, SWORK( PTSX ), N, X, LDX, INFO );

      // Compute R = B - AX (R is WORK).

      dlacpy('All', N, NRHS, B, LDB, WORK, N );

      dgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER=0 and return.

      for (I = 1; I <= NRHS; I++) {
         XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
         RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
         if (RNRM > XNRM*CTE) GO TO 10;
      }

      // If we are here, the NRHS normwise backward errors satisfy the
      // stopping criterion. We are good to exit.

      ITER = 0
      RETURN

      } // 10

      for (IITER = 1; IITER <= ITERMAX; IITER++) { // 30

         // Convert R (in WORK) from double precision to single precision
         // and store the result in SX.

         dlag2s(N, NRHS, WORK, N, SWORK( PTSX ), N, INFO );

         if ( INFO != 0 ) {
            ITER = -2
            GO TO 40
         }

         // Solve the system SA*SX = SR.

         sgetrs('No transpose', N, NRHS, SWORK( PTSA ), N, IPIV, SWORK( PTSX ), N, INFO );

         // Convert SX back to double precision and update the current
         // iterate.

         slag2d(N, NRHS, SWORK( PTSX ), N, WORK, N, INFO );

         for (I = 1; I <= NRHS; I++) {
            daxpy(N, ONE, WORK( 1, I ), 1, X( 1, I ), 1 );
         }

         // Compute R = B - AX (R is WORK).

         dlacpy('All', N, NRHS, B, LDB, WORK, N );

         dgemm('No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

         // Check whether the NRHS normwise backward errors satisfy the
         // stopping criterion. If yes, set ITER=IITER>0 and return.

         for (I = 1; I <= NRHS; I++) {
            XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
            RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
            if (RNRM > XNRM*CTE) GO TO 20;
         }

         // If we are here, the NRHS normwise backward errors satisfy the
         // stopping criterion, we are good to exit.

         ITER = IITER

         RETURN

         } // 20

      } // 30

      // If we are at this place of the code, this is because we have
      // performed ITER=ITERMAX iterations and never satisfied the
      // stopping criterion, set up the ITER flag accordingly and follow up
      // on double precision routine.

      ITER = -ITERMAX - 1

      } // 40

      // Single-precision iterative refinement failed to converge to a
      // satisfactory solution, so we resort to double precision.

      dgetrf(N, N, A, LDA, IPIV, INFO );

      if (INFO != 0) RETURN;

      dlacpy('All', N, NRHS, B, LDB, X, LDX );
      dgetrs('No transpose', N, NRHS, A, LDA, IPIV, X, LDX, INFO );

      RETURN

      // End of DSGESV

      }
