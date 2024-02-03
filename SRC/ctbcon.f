      SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, KD, LDAB, N;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      REAL               AINVNM, ANORM, SCALE, SMLNUM, XNORM
      COMPLEX            ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               CLANTB, SLAMCH
      // EXTERNAL LSAME, ICAMAX, CLANTB, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATBS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      ONENRM = NORM == '1' || LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )

      if ( !ONENRM && !LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( KD < 0 ) {
         INFO = -5
      } else if ( LDAB < KD+1 ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('CTBCON', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         RCOND = ONE
         RETURN
      }

      RCOND = ZERO
      SMLNUM = SLAMCH( 'Safe minimum' )*REAL( MAX( N, 1 ) )

      // Compute the 1-norm of the triangular matrix A or A**H.

      ANORM = CLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, RWORK )

      // Continue only if ANORM > 0.

      if ( ANORM > ZERO ) {

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
         clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == KASE1 ) {

               // Multiply by inv(A).

               clatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO );
            } else {

               // Multiply by inv(A**H).

               clatbs(UPLO, 'Conjugate transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, RWORK, INFO );
            }
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE != ONE ) {
               IX = ICAMAX( N, WORK, 1 )
               XNORM = CABS1( WORK( IX ) )
               if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 20;
               csrscl(N, SCALE, WORK, 1 );
            }
            GO TO 10
         }

         // Compute the estimate of the reciprocal condition number.

         if (AINVNM != ZERO) RCOND = ( ONE / ANORM ) / AINVNM;
      }

      } // 20
      RETURN

      // End of CTBCON

      }
