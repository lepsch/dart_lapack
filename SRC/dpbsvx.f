      SUBROUTINE DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), FERR( * ), S( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU, UPPER;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSB;
      // EXTERNAL LSAME, DLAMCH, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLACPY, DLAQSB, DPBCON, DPBEQU, DPBRFS, DPBTRF, DPBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      UPPER = LSAME( UPLO, 'U' )
      if ( NOFACT .OR. EQUIL ) {
         EQUED = 'N'
         RCEQU = false;
      } else {
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      }

      // Test the input parameters.

      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KD.LT.0 ) {
         INFO = -4
      } else if ( NRHS.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -7
      } else if ( LDAFB.LT.KD+1 ) {
         INFO = -9
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -10
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM
            SMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
            } // 10
            if ( SMIN.LE.ZERO ) {
               INFO = -11
            } else if ( N.GT.0 ) {
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            }
         }
         if ( INFO.EQ.0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -13
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -15
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DPBSVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         dpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

            // Equilibrate the matrix.

            dlaqsb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' )
         }
      }

      // Scale the right-hand side.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 30
            for (I = 1; I <= N; I++) { // 20
               B( I, J ) = S( I )*B( I, J )
            } // 20
         } // 30
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the Cholesky factorization A = U**T *U or A = L*L**T.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 40
               J1 = MAX( J-KD, 1 )
               dcopy(J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 );
            } // 40
         } else {
            for (J = 1; J <= N; J++) { // 50
               J2 = MIN( J+KD, N )
               dcopy(J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 );
            } // 50
         }

         dpbtrf(UPLO, N, KD, AFB, LDAFB, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = DLANSB( '1', UPLO, N, KD, AB, LDAB, WORK )

      // Compute the reciprocal of the condition number of A.

      dpbcon(UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution matrix X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dpbtrs(UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      dpbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 70
            for (I = 1; I <= N; I++) { // 60
               X( I, J ) = S( I )*X( I, J )
            } // 60
         } // 70
         for (J = 1; J <= NRHS; J++) { // 80
            FERR( J ) = FERR( J ) / SCOND
         } // 80
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of DPBSVX

      }
