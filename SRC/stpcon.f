      SUBROUTINE STPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, N;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AP( * ), WORK( * )
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
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH, SLANTP
      // EXTERNAL LSAME, ISAMAX, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATPS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
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
      }
      if ( INFO != 0 ) {
         xerbla('STPCON', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         RCOND = ONE
         RETURN
      }

      RCOND = ZERO
      SMLNUM = SLAMCH( 'Safe minimum' )*REAL( MAX( 1, N ) )

      // Compute the norm of the triangular matrix A.

      ANORM = SLANTP( NORM, UPLO, DIAG, N, AP, WORK )

      // Continue only if ANORM > 0.

      if ( ANORM > ZERO ) {

         // Estimate the norm of the inverse of A.

         AINVNM = ZERO
         NORMIN = 'N'
         if ( ONENRM ) {
            KASE1 = 1
         } else {
            KASE1 = 2
         }
         KASE = 0
         } // 10
         slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == KASE1 ) {

               // Multiply by inv(A).

               slatps(UPLO, 'No transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO );
            } else {

               // Multiply by inv(A**T).

               slatps(UPLO, 'Transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO );
            }
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE != ONE ) {
               IX = ISAMAX( N, WORK, 1 )
               XNORM = ABS( WORK( IX ) )
               if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 20;
               srscl(N, SCALE, WORK, 1 );
            }
            GO TO 10
         }

         // Compute the estimate of the reciprocal condition number.

         if (AINVNM != ZERO) RCOND = ( ONE / ANORM ) / AINVNM;
      }

      } // 20
      RETURN

      // End of STPCON

      }
