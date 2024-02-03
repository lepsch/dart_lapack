      void dtbcon(NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, KD, LDAB, N;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, ANORM, SCALE, SMLNUM, XNORM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH, DLANTB;
      // EXTERNAL LSAME, IDAMAX, DLAMCH, DLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATBS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      ONENRM = NORM == '1' || LSAME( NORM, 'O' );
      NOUNIT = LSAME( DIAG, 'N' );

      if ( !ONENRM && !LSAME( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( KD < 0 ) {
         INFO = -5;
      } else if ( LDAB < KD+1 ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DTBCON', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         RCOND = ONE;
         return;
      }

      RCOND = ZERO;
      SMLNUM = DLAMCH( 'Safe minimum' )*DBLE( max( 1, N ) );

      // Compute the norm of the triangular matrix A.

      ANORM = DLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, WORK );

      // Continue only if ANORM > 0.

      if ( ANORM > ZERO ) {

         // Estimate the norm of the inverse of A.

         AINVNM = ZERO;
         NORMIN = 'N';
         if ( ONENRM ) {
            KASE1 = 1;
         } else {
            KASE1 = 2;
         }
         KASE = 0;
         } // 10
         dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == KASE1 ) {

               // Multiply by inv(A).

               dlatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
            } else {

               // Multiply by inv(A**T).

               dlatbs(UPLO, 'Transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
            }
            NORMIN = 'Y';

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE != ONE ) {
               IX = IDAMAX( N, WORK, 1 );
               XNORM = ABS( WORK( IX ) );
               if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 20;
               drscl(N, SCALE, WORK, 1 );
            }
            GO TO 10;
         }

         // Compute the estimate of the reciprocal condition number.

         if (AINVNM != ZERO) RCOND = ( ONE / ANORM ) / AINVNM;
      }

      } // 20
      return;

      // End of DTBCON

      }
