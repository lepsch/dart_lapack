      SUBROUTINE ZCGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, SWORK, RWORK, ITER, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, ITER, LDA, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             RWORK( * );
      COMPLEX            SWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ), X( LDX, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      bool               DOITREF;
      PARAMETER          ( DOITREF = .TRUE. )
*
      int                ITERMAX;
      PARAMETER          ( ITERMAX = 30 )
*
      double             BWDMAX;
      PARAMETER          ( BWDMAX = 1.0E+00 )
*
      COMPLEX*16         NEGONE, ONE
      PARAMETER          ( NEGONE = ( -1.0D+00, 0.0D+00 ), ONE = ( 1.0D+00, 0.0D+00 ) )
*
      // .. Local Scalars ..
      int                I, IITER, PTSA, PTSX;
      double             ANRM, CTE, EPS, RNRM, XNRM;
      COMPLEX*16         ZDUM
*
      // .. External Subroutines ..
      // EXTERNAL CGETRS, CGETRF, CLAG2Z, XERBLA, ZAXPY, ZGEMM, ZLACPY, ZLAG2C, ZGETRF, ZGETRS
      // ..
      // .. External Functions ..
      int                IZAMAX;
      double             DLAMCH, ZLANGE;
      // EXTERNAL IZAMAX, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      ITER = 0
*
      // Test the input parameters.
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZCGESV', -INFO )
         RETURN
      END IF
*
      // Quick return if (N.EQ.0).
*
      IF( N.EQ.0 ) RETURN
*
      // Skip single precision iterative refinement if a priori slower
     t // han double precision factorization.
*
      IF( .NOT.DOITREF ) THEN
         ITER = -1
         GO TO 40
      END IF
*
      // Compute some constants.
*
      ANRM = ZLANGE( 'I', N, N, A, LDA, RWORK )
      EPS = DLAMCH( 'Epsilon' )
      CTE = ANRM*EPS*SQRT( DBLE( N ) )*BWDMAX
*
      // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.
*
      PTSA = 1
      PTSX = PTSA + N*N
*
      // Convert B from double precision to single precision and store the
      // result in SX.
*
      CALL ZLAG2C( N, NRHS, B, LDB, SWORK( PTSX ), N, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -2
         GO TO 40
      END IF
*
      // Convert A from double precision to single precision and store the
      // result in SA.
*
      CALL ZLAG2C( N, N, A, LDA, SWORK( PTSA ), N, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -2
         GO TO 40
      END IF
*
      // Compute the LU factorization of SA.
*
      CALL CGETRF( N, N, SWORK( PTSA ), N, IPIV, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -3
         GO TO 40
      END IF
*
      // Solve the system SA*SX = SB.
*
      CALL CGETRS( 'No transpose', N, NRHS, SWORK( PTSA ), N, IPIV, SWORK( PTSX ), N, INFO )
*
      // Convert SX back to double precision
*
      CALL CLAG2Z( N, NRHS, SWORK( PTSX ), N, X, LDX, INFO )
*
      // Compute R = B - AX (R is WORK).
*
      CALL ZLACPY( 'All', N, NRHS, B, LDB, WORK, N )
*
      CALL ZGEMM( 'No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX, ONE, WORK, N )
*
      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER=0 and return.
*
      DO I = 1, NRHS
         XNRM = CABS1( X( IZAMAX( N, X( 1, I ), 1 ), I ) )
         RNRM = CABS1( WORK( IZAMAX( N, WORK( 1, I ), 1 ), I ) )
         IF( RNRM.GT.XNRM*CTE ) GO TO 10
      END DO
*
      // If we are here, the NRHS normwise backward errors satisfy the
      // stopping criterion. We are good to exit.
*
      ITER = 0
      RETURN
*
   10 CONTINUE
*
      DO 30 IITER = 1, ITERMAX
*
         // Convert R (in WORK) from double precision to single precision
         // and store the result in SX.
*
         CALL ZLAG2C( N, NRHS, WORK, N, SWORK( PTSX ), N, INFO )
*
         IF( INFO.NE.0 ) THEN
            ITER = -2
            GO TO 40
         END IF
*
         // Solve the system SA*SX = SR.
*
         CALL CGETRS( 'No transpose', N, NRHS, SWORK( PTSA ), N, IPIV, SWORK( PTSX ), N, INFO )
*
         // Convert SX back to double precision and update the current
         // iterate.
*
         CALL CLAG2Z( N, NRHS, SWORK( PTSX ), N, WORK, N, INFO )
*
         DO I = 1, NRHS
            CALL ZAXPY( N, ONE, WORK( 1, I ), 1, X( 1, I ), 1 )
         END DO
*
         // Compute R = B - AX (R is WORK).
*
         CALL ZLACPY( 'All', N, NRHS, B, LDB, WORK, N )
*
         CALL ZGEMM( 'No Transpose', 'No Transpose', N, NRHS, N, NEGONE, A, LDA, X, LDX, ONE, WORK, N )
*
         // Check whether the NRHS normwise backward errors satisfy the
         // stopping criterion. If yes, set ITER=IITER>0 and return.
*
         DO I = 1, NRHS
            XNRM = CABS1( X( IZAMAX( N, X( 1, I ), 1 ), I ) )
            RNRM = CABS1( WORK( IZAMAX( N, WORK( 1, I ), 1 ), I ) )
            IF( RNRM.GT.XNRM*CTE ) GO TO 20
         END DO
*
         // If we are here, the NRHS normwise backward errors satisfy the
         // stopping criterion, we are good to exit.
*
         ITER = IITER
*
         RETURN
*
   20    CONTINUE
*
   30 CONTINUE
*
      // If we are at this place of the code, this is because we have
      // performed ITER=ITERMAX iterations and never satisfied the stopping
      // criterion, set up the ITER flag accordingly and follow up on double
      // precision routine.
*
      ITER = -ITERMAX - 1
*
   40 CONTINUE
*
      // Single-precision iterative refinement failed to converge to a
      // satisfactory solution, so we resort to double precision.
*
      CALL ZGETRF( N, N, A, LDA, IPIV, INFO )
*
      IF( INFO.NE.0 ) RETURN
*
      CALL ZLACPY( 'All', N, NRHS, B, LDB, X, LDX )
      CALL ZGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, X, LDX, INFO )
*
      RETURN
*
      // End of ZCGESV
*
      END
