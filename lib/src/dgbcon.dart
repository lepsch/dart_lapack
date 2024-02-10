import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgbcon(NORM, N, KL, KU, final Matrix<double> AB, final int LDAB, final Array<int> IPIV, ANORM, RCOND, final Array<double> _WORK, final Array<int> IWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                INFO, KL, KU, LDAB, N;
      double             ANORM, RCOND;
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               LNOTI, ONENRM;
      String             NORMIN;
      int                IX, J, JP, KASE, KASE1, KD, LM;
      double             AINVNM, SCALE, SMLNUM, T;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DDOT, DLAMCH;
      // EXTERNAL lsame, idamax, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLACN2, DLATBS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN

      // Test the input parameters.

      INFO = 0;
      ONENRM = NORM == '1' || lsame( NORM, 'O' );
      if ( !ONENRM && !lsame( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 ) {
         INFO = -3;
      } else if ( KU < 0 ) {
         INFO = -4;
      } else if ( LDAB < 2*KL+KU+1 ) {
         INFO = -6;
      } else if ( ANORM < ZERO ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DGBCON', -INFO );
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

      SMLNUM = dlamch( 'Safe minimum' );

      // Estimate the norm of inv(A).

      AINVNM = ZERO;
      NORMIN = 'N';
      if ( ONENRM ) {
         KASE1 = 1;
      } else {
         KASE1 = 2;
      }
      KD = KL + KU + 1;
      LNOTI = KL > 0;
      KASE = 0;
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            if ( LNOTI ) {
               for (J = 1; J <= N - 1; J++) { // 20
                  LM = min( KL, N-J );
                  JP = IPIV( J );
                  T = WORK( JP );
                  if ( JP != J ) {
                     WORK[JP] = WORK( J );
                     WORK[J] = T;
                  }
                  daxpy(LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 );
               } // 20
            }

            // Multiply by inv(U).

            dlatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(U**T).

            dlatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );

            // Multiply by inv(L**T).

            if ( LNOTI ) {
               for (J = N - 1; J >= 1; J--) { // 30
                  LM = min( KL, N-J );
                  WORK[J] = WORK( J ) - ddot( LM, AB( KD+1, J ), 1, WORK( J+1 ), 1 );
                  JP = IPIV( J );
                  if ( JP != J ) {
                     T = WORK( JP );
                     WORK[JP] = WORK( J );
                     WORK[J] = T;
                  }
               } // 30
            }
         }

         // Divide X by 1/SCALE if doing so will not cause overflow.

         NORMIN = 'Y';
         if ( SCALE != ONE ) {
            IX = idamax( N, WORK, 1 );
            if( SCALE < ( WORK( IX ) ).abs()*SMLNUM || SCALE == ZERO ) GO TO 40;
            drscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 40
      }
