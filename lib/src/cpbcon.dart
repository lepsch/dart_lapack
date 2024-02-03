      void cpbcon(UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      REAL               ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            AB( LDAB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM;
      COMPLEX            ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATBS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ( REAL( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDAB < KD+1 ) {
         INFO = -5;
      } else if ( ANORM < ZERO ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CPBCON', -INFO );
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

            clatbs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y';

            // Multiply by inv(U).

            clatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO );
         } else {

            // Multiply by inv(L).

            clatbs('Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK, INFO );
            NORMIN = 'Y';

            // Multiply by inv(L**H).

            clatbs('Lower', 'Conjugate transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK, INFO );
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

      return;
      }
