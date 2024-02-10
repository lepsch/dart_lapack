      void cppcon(UPLO, N, AP, ANORM, RCOND, WORK, final Array<double> RWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      double               ANORM, RCOND;
      double               RWORK( * );
      Complex            AP( * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      double               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM;
      Complex            ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATPS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ANORM < ZERO ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CPPCON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      }

      SMLNUM = SLAMCH( 'Safe minimum' );

      // Estimate the 1-norm of the inverse.

      KASE = 0;
      NORMIN = 'N';
      } // 10
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**H).

            clatps('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y';

            // Multiply by inv(U).

            clatps('Upper', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEU, RWORK, INFO );
         } else {

            // Multiply by inv(L).

            clatps('Lower', 'No transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y';

            // Multiply by inv(L**H).

            clatps('Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, AP, WORK, SCALEU, RWORK, INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU;
         if ( SCALE != ONE ) {
            IX = ICAMAX( N, WORK, 1 );
            if( SCALE < CABS1( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20;
            csrscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 20
      }
