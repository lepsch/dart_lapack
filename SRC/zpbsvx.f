      SUBROUTINE ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RWORK( * ), S( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
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
      double             DLAMCH, ZLANHB;
      // EXTERNAL LSAME, DLAMCH, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAQHB, ZPBCON, ZPBEQU, ZPBRFS, ZPBTRF, ZPBTRS
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
         RCEQU = .FALSE.
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
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
   10       CONTINUE
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
         xerbla('ZPBSVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         zpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

            // Equilibrate the matrix.

            zlaqhb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' )
         }
      }

      // Scale the right-hand side.

      if ( RCEQU ) {
         DO 30 J = 1, NRHS
            DO 20 I = 1, N
               B( I, J ) = S( I )*B( I, J )
   20       CONTINUE
   30    CONTINUE
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

         if ( UPPER ) {
            DO 40 J = 1, N
               J1 = MAX( J-KD, 1 )
               zcopy(J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 );
   40       CONTINUE
         } else {
            DO 50 J = 1, N
               J2 = MIN( J+KD, N )
               zcopy(J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 );
   50       CONTINUE
         }

         zpbtrf(UPLO, N, KD, AFB, LDAFB, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANHB( '1', UPLO, N, KD, AB, LDAB, RWORK )

      // Compute the reciprocal of the condition number of A.

      zpbcon(UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, RWORK, INFO );

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zpbtrs(UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zpbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         DO 70 J = 1, NRHS
            DO 60 I = 1, N
               X( I, J ) = S( I )*X( I, J )
   60       CONTINUE
   70    CONTINUE
         DO 80 J = 1, NRHS
            FERR( J ) = FERR( J ) / SCOND
   80    CONTINUE
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of ZPBSVX

      }
