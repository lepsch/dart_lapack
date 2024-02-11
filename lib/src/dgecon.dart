import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgecon(final int NORM, final int N, final Matrix<double> A, final int LDA, final int ANORM, final int RCOND, final Array<double> _WORK, final Array<int> IWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      int                IWORK( * );
      double             A( LDA, * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      //- int                idamax;
      //- double             DLAMCH;
      // EXTERNAL lsame, idamax, DLAMCH, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATRS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      HUGEVAL = dlamch( 'Overflow' );

      // Test the input parameters.

      INFO = 0;
      ONENRM = NORM == '1' || lsame( NORM, 'O' );
      if ( !ONENRM && !lsame( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('DGECON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      } else if ( disnan( ANORM ) ) {
         RCOND = ANORM;
         INFO = -5;
         return;
      } else if ( ANORM > HUGEVAL ) {
         INFO = -5;
         return;
      }

      SMLNUM = dlamch( 'Safe minimum' );

      // Estimate the norm of inv(A).

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

            // Multiply by inv(L).

            dlatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );

            // Multiply by inv(U).

            dlatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );
         } else {

            // Multiply by inv(U**T).

            dlatrs('Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );

            // Multiply by inv(L**T).

            dlatrs('Lower', 'Transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );
         }

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU;
         NORMIN = 'Y';
         if ( SCALE != ONE ) {
            IX = idamax( N, WORK, 1 );
            if( SCALE < ( WORK( IX ) ).abs()*SMLNUM || SCALE == ZERO ) GO TO 20;
            drscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if ( AINVNM != ZERO ) {
         RCOND = ( ONE / AINVNM ) / ANORM;
      } else {
         INFO = 1;
         return;
      }

      // Check for NaNs and Infs

      if( disnan( RCOND ) || RCOND > HUGEVAL ) INFO = 1;

      } // 20
      }
