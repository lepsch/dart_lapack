import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtpcon(NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, NORM, UPLO;
      int                INFO, N;
      double             RCOND;
      int                IWORK( * );
      double             AP( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               NOUNIT, ONENRM, UPPER;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, ANORM, SCALE, SMLNUM, XNORM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DLAMCH, DLANTP;
      // EXTERNAL lsame, idamax, DLAMCH, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATPS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      ONENRM = NORM == '1' || lsame( NORM, 'O' );
      NOUNIT = lsame( DIAG, 'N' );

      if ( !ONENRM && !lsame( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DTPCON', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         RCOND = ONE;
         return;
      }

      RCOND = ZERO;
      SMLNUM = dlamch( 'Safe minimum' )*(max( 1, N )).toDouble();

      // Compute the norm of the triangular matrix A.

      ANORM = DLANTP( NORM, UPLO, DIAG, N, AP, WORK );

      // Continue only if ANORM > 0.

      if ( ANORM > ZERO ) {

         // Estimate the norm of the inverse of A.

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

               // Multiply by inv(A).

               dlatps(UPLO, 'No transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO );
            } else {

               // Multiply by inv(A**T).

               dlatps(UPLO, 'Transpose', DIAG, NORMIN, N, AP, WORK, SCALE, WORK( 2*N+1 ), INFO );
            }
            NORMIN = 'Y';

            // Multiply by 1/SCALE if doing so will not cause overflow.

            if ( SCALE != ONE ) {
               IX = idamax( N, WORK, 1 );
               XNORM = ( WORK( IX ) ).abs();
               if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 20;
               drscl(N, SCALE, WORK, 1 );
            }
            GO TO 10;
         }

         // Compute the estimate of the reciprocal condition number.

         if (AINVNM != ZERO) RCOND = ( ONE / ANORM ) / AINVNM;
      }

      } // 20
      }
