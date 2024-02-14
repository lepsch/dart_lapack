      void spbcon(final int UPLO, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final int ANORM, final int RCOND, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final AB = AB_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, N;
      double               ANORM, RCOND;
      int                IWORK( * );
      double               AB( LDAB, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      double               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ISAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATBS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
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
         xerbla('SPBCON', -INFO );
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
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**T).

            slatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y';

            // Multiply by inv(U).

            slatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(L).

            slatbs('Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y';

            // Multiply by inv(L**T).

            slatbs('Lower', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU;
         if ( SCALE != ONE ) {
            IX = ISAMAX( N, WORK, 1 );
            if( SCALE < ( WORK( IX ) ).abs()*SMLNUM || SCALE == ZERO ) GO TO 20;
            srscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 20

      }
