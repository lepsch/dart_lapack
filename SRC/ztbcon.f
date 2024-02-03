      SUBROUTINE ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, KD, LDAB, N;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         AB( LDAB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, ANORM, SCALE, SMLNUM, XNORM;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, ZLANTB;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, ZLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDRSCL, ZLACN2, ZLATBS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      ONENRM = NORM == '1' .OR. LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )

      if ( .NOT.ONENRM && .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT && .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('ZTBCON', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         RCOND = ONE
         RETURN
      }

      RCOND = ZERO
      SMLNUM = DLAMCH( 'Safe minimum' )*DBLE( MAX( N, 1 ) )

      // Compute the 1-norm of the triangular matrix A or A**H.

      ANORM = ZLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK )

      // Continue only if ANORM > 0.

      if ( ANORM.GT.ZERO ) {

         // Estimate the 1-norm of the inverse of A.

         AINVNM = ZERO
         NORMIN = 'N'
         if ( ONENRM ) {
            KASE1 = 1
         } else {
            KASE1 = 2
         }
         KASE = 0
         } // 10
         zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == KASE1 ) {

               // Multiply by inv(A).

               zlatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO );
            } else {

               // Multiply by inv(A**H).

               zlatbs(UPLO, 'Conjugate transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO );
            }
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE != ONE ) {
               IX = IZAMAX( N, WORK, 1 )
               XNORM = CABS1( WORK( IX ) )
               if (SCALE.LT.XNORM*SMLNUM .OR. SCALE == ZERO) GO TO 20;
               zdrscl(N, SCALE, WORK, 1 );
            }
            GO TO 10
         }

         // Compute the estimate of the reciprocal condition number.

         if (AINVNM != ZERO) RCOND = ( ONE / ANORM ) / AINVNM;
      }

      } // 20
      RETURN

      // End of ZTBCON

      }
