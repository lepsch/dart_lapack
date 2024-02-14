import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtbcon(final int NORM, final int UPLO, final int DIAG, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final int RCOND, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final AB = AB_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, NORM, UPLO;
      int                INFO, KD, LDAB, N;
      double             RCOND;
      int                IWORK( * );
      double             AB( LDAB, * ), WORK( * );
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
      //- double             DLAMCH, DLANTB;
      // EXTERNAL lsame, idamax, DLAMCH, DLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATBS, DRSCL, XERBLA
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
      } else if ( KD < 0 ) {
         INFO = -5;
      } else if ( LDAB < KD+1 ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DTBCON', -INFO );
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

      ANORM = DLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, WORK );

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

               dlatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
            } else {

               // Multiply by inv(A**T).

               dlatbs(UPLO, 'Transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
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
