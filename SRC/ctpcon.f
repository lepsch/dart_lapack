      SUBROUTINE CTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                INFO, N;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AP( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
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
      REAL               CLANTP, SLAMCH
      // EXTERNAL LSAME, ICAMAX, CLANTP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATPS, CSRSCL, XERBLA
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
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )

      if ( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('CTPCON', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      }

      RCOND = ZERO
      SMLNUM = SLAMCH( 'Safe minimum' )*REAL( MAX( 1, N ) )

      // Compute the norm of the triangular matrix A.

      ANORM = CLANTP( NORM, UPLO, DIAG, N, AP, RWORK )

      // Continue only if ANORM > 0.

      if ( ANORM.GT.ZERO ) {

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
         clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.KASE1 ) {

               // Multiply by inv(A).

               clatps(UPLO, 'No transpose', DIAG, NORMIN, N, AP, WORK, SCALE, RWORK, INFO );
            } else {

               // Multiply by inv(A**H).

               clatps(UPLO, 'Conjugate transpose', DIAG, NORMIN, N, AP, WORK, SCALE, RWORK, INFO );
            }
            NORMIN = 'Y'

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE.NE.ONE ) {
               IX = ICAMAX( N, WORK, 1 )
               XNORM = CABS1( WORK( IX ) )
               IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
               csrscl(N, SCALE, WORK, 1 );
            }
            GO TO 10
         }

         // Compute the estimate of the reciprocal condition number.

         IF( AINVNM.NE.ZERO ) RCOND = ( ONE / ANORM ) / AINVNM
      }

      } // 20
      RETURN

      // End of CTPCON

      }
