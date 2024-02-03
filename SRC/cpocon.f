      SUBROUTINE CPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      COMPLEX            ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATRS, CSRSCL, XERBLA
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
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( ANORM.LT.ZERO ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('CPOCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM.EQ.ZERO ) {
         RETURN
      }

      SMLNUM = SLAMCH( 'Safe minimum' )

      // Estimate the 1-norm of inv(A).

      KASE = 0
      NORMIN = 'N'
      } // 10
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**H).

            clatrs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y'

            // Multiply by inv(U).

            clatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEU, RWORK, INFO );
         } else {

            // Multiply by inv(L).

            clatrs('Lower', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y'

            // Multiply by inv(L**H).

            clatrs('Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEU, RWORK, INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU
         if ( SCALE.NE.ONE ) {
            IX = ICAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            csrscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      } // 20
      RETURN

      // End of CPOCON

      }
